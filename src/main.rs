#![allow(unused_parens)]

extern crate crossbeam;

use crossbeam::queue::PopError;
use crossbeam::queue::SegQueue;
use std::collections::BTreeSet;
use std::collections::HashMap;

mod flags;

use flags::Flags;

fn egcd(a: isize, b: isize) -> (isize, isize, isize) {
    let mut a = a.abs();
    let mut b = b.abs();

    let mut sa = 1;
    let mut ta = 0;

    let mut sb = 0;
    let mut tb = 1;

    while a > 0 {
        let q = b / a;

        sb -= sa * q;
        tb -= ta * q;
        b -= a * q;

        std::mem::swap(&mut sa, &mut sb);
        std::mem::swap(&mut ta, &mut tb);
        std::mem::swap(&mut a, &mut b);
    }

    (b, sb, tb)
}

fn main() {
    let n = 12;

    let threads = 8;
    let workunit_bits = 6;

    let mut lattices = Vec::new();
    for w in 1..=n {
        if n % w != 0 {
            continue;
        }

        let h = n / w;

        for sx in 0..w {
            if 2 * sx > w {
                // don't bother with shifts over half, instead we'll find whatever flipped
                continue;
            }

            let (d, s, t) = egcd(sx, w);
            // d = s * sx + t * w

            // (w, 0), (sx, h)
            // (w, 0), (sx, h), (sx * w / d, h * w / d), (s * sx + t * w, s * h + t * 0)
            // (w, 0), (sx, h), (sx/d * w, h * w / d), (d, s * h)
            // (w, 0), (sx, h), (0, h * w / d), (d, s * h)

            // the order stepping vertically, what w would be if we transposed
            let wt = h * w / d;
            // the rest...
            let ht = d;
            // what sx would be if we transposed
            let sxt = s * t;
            // grumble grumble stupid mod
            let sxt = ((sxt % wt) + wt) % wt;
            // if smaller flipped...
            let sxt = sxt.min(wt - sxt);

            if (-wt, ht, sxt) < (-w, ht, sx) {
                continue;
            }

            lattices.push((w, h, sx));
        }
    }
    //println!("Lattices {:?}", lattices);

    for lattice in lattices {
        let (w, h, _sx) = lattice;

        let flags = Flags::new(1 << n);
        let flags = &flags;

        let workunits: Vec<_> = (0..(1 << workunit_bits)).collect();
        let mut results: Vec<_> = workunits.iter().map(|_| Vec::new()).collect();

        {
            let q = SegQueue::new();
            for tuple in workunits.into_iter().zip(results.iter_mut()) {
                q.push(tuple);
            }

            crossbeam::scope(|sc| {
                for _ in 0..threads {
                    sc.spawn(|_| {
                        loop {
                            let (workunit, results) = match q.pop() {
                                Result::Ok(tuple) => tuple,
                                Result::Err(PopError) => {
                                    return;
                                }
                            };

                            let suffix_bits = n - workunit_bits;
                            for suffix in 0..(1 << suffix_bits) {
                                let s0 = (workunit << suffix_bits) | suffix;
                                search(lattice, flags, s0, results);
                            }
                        }
                    });
                }
            }).unwrap();
        }

        let results = results.into_iter().flatten();
        let results = results.filter(|&(s, _)| {
            for dx in 0..w {
                for dy in 0..h {
                    if dx == 0 && dy == 0 {
                        continue;
                    }

                    let mut s2 = 0;
                    for x in 0..w {
                        for y in 0..h {
                            let idx = y * w + x;

                            let (x2, y2) = canon_2d(lattice, x + dx, y + dy);
                            let idx2 = y2 * w + x2;

                            if (s >> idx) & 1 == 1 {
                                s2 |= (1 << idx2);
                            }
                        }
                    }

                    if s == s2 {
                        //eprintln!("Drop lattice {:?} result {} shift ({}, {}), dx, dy", lattice, s, dx, dy);
                        return false;
                    }
                }
            }

            true
        });
        let results: BTreeSet<_> = results.collect();
        eprintln!("Lattice {:?} => {} results", lattice, results.len());
        for result in results {
            eprintln!("   {:?}", result);
        }
    }

    // TODO: more analysis...

    // analyze connected components in space and time
    //
    // this is directed graph of (x, y, t) with edges labelled with lattice shift (lattice shift is
    // number of x wraps, number of y wraps, number of t wraps)
    //
    // all our cells (in t = 0) had better be part of same connected component or we discard
    // (should find connected components separately).
    //
    // then find cycles whose edge label sums tell us generators for "lattice of lattice shifts"
    //
    // These generators can be recast as actual x/y/t distances to copies of same cell that it's
    // connected to.  Either of these lattices are the same for further purposes although latter is
    // more human-intelligible I think.

    // rank of intersection of this connection lattice with t = 0 tells us...
    //
    // zero rank: Oscillator or glider, probably discard since we don't expect any interesting
    // results.  Could analyze as oscillator/glider to give period and shift.
    //
    // one rank: Wick of some sort.  Presumably all interesting although overpop-only connection
    // may mean a lot of boring stuff here.  Projection of connection lattice into (x, y) plane
    // (rather than intersection with t = 0) tells us stuff here.  If it is also one rank then
    // we've got an oscillator wick and if it's two rank we have a (sideways) moving wick.
    //
    // two rank: Real agar.
}

fn search(lattice: (isize, isize, isize), flags: &Flags, s0: u64, results: &mut Vec<(u64, usize)>) {
    let mut prev_vec = Vec::new();
    let mut prev_map = HashMap::new();
    let mut s = s0;

    loop {
        if flags.get(s) {
            break;
        }

        if let Some(&idx) = prev_map.get(&s) {
            let &sc = (&prev_vec[idx..]).iter().min().unwrap();
            results.push((sc, prev_vec.len() - idx));
            break;
        }

        let t = prev_vec.len();
        prev_vec.push(s);
        prev_map.insert(s, t);
        s = tick(lattice, s);
    }

    for sp in prev_vec {
        flags.set(sp);
    }
}

fn tick(lattice: (isize, isize, isize), s0: u64) -> u64 {
    let (w, h, _sx) = lattice;
    let mut s1 = 0;
    for x in 0..w {
        for y in 0..h {
            let idx = y * w + x;
            let mut s = 0;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    let (x2, y2) = canon_2d(lattice, x + dx, y + dy);
                    let idx2 = y2 * w + x2;
                    s += ((s0 >> idx2) & 1);
                }
            }
            let living = match ((s0 >> idx) & 1 == 1) {
                true => (3 <= s && s <= 4),
                false => (s == 3),
            };
            if living {
                s1 |= (1 << idx);
            }
        }
    }
    s1
}

fn canon_2d(lattice: (isize, isize, isize), x: isize, y: isize) -> (isize, isize) {
    let (w, h, sx) = lattice;
    let mut x = x;
    let mut y = y;

    while y < 0 {
        x += sx;
        y += h;
    }
    while y >= h {
        x -= sx;
        y -= h;
    }
    while x < 0 {
        x += w;
    }
    while x >= w {
        x -= w;
    }
    (x, y)
}
