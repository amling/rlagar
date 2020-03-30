#![allow(unused_parens)]

extern crate crossbeam;

use crossbeam::queue::PopError;
use crossbeam::queue::SegQueue;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;

mod algo;
mod flags;
mod lattice;

use flags::Flags;
use lattice::LatticeCanonicalizable;

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

            let (d, _, (s, t)) = lattice::egcd(sx, w, (1, 0), (0, 1));
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
        let (w, h, sx) = lattice;

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

        // drop results that have an additional symmetry since they'll be found in smaller
        // geometries
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

                            let (x2, y2, _, _) = canon_2d(lattice, x + dx, y + dy);
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

        // TODO: Drop results which have a smaller shift (in any generation).  We're getting tons
        // of duplicates that are just shifts of each other.

        let results: BTreeSet<_> = results.collect();
        eprintln!("Lattice {:?} => {} results", lattice, results.len());
        for result in results {
            //eprintln!("   {:?}", result);

            let (s, period) = result;
            let mut links = HashMap::new();
            let mut s1 = s;
            for t in 0..period {
                for ((x1, y1, t1), links2) in compute_links(lattice, s1).into_iter() {
                    let t1 = t + t1;
                    let (t1, lt1) = canon_1d(period, t1);
                    for ((x2, y2, t2), (lx, ly)) in links2.into_iter() {
                        let t2 = t + t2;
                        let (t2, lt2) = canon_1d(period, t2);
                        links.entry((x1, y1, t1)).or_insert_with(|| HashSet::new()).insert(((x2, y2, t2), (lx, ly, lt2 - lt1)));
                    }
                }
                s1 = tick(lattice, s1);
            }

            // links is the directed graph of (x, y, t) with edges labelled with lattice shift
            // (lattice shift is number of x wraps, number of y wraps, number of t wraps).

            if let Some(ls) = compute_lattice_links(&links) {
                // These "lattice links" tell us the generators for the "lattice of lattice shifts"
                // that connect the same cell in each fundamntal period (in the 3D lattice of space
                // and time).  These are so far given in terms of wrap counts, not actualy x/y
                // coordinates.

                // These generators can be recast as actual x/y/t distances to copies of same cell
                // that it's connected to.  Either of these lattices are the same for further
                // purposes (rank analysis) although latter is more human-intelligible I think.
                let ls: Vec<_> = ls.into_iter().map(|(lx, ly, lt)| (lx * w + ly * sx, ly * h, lt * period)).collect();

                // now reform as ((((), x), y), t) for lattice canonicalizer
                let fl: Vec<_> = ls.iter().map(|&(x, y, t)| ((((), x), y), t)).collect();
                let fl = LatticeCanonicalizable::canonicalize(fl);

                let (fl_vt, fl) = fl;
                let fl_vt = fl_vt.unwrap();
                let fl_vt = (((fl_vt.0).0).1, (fl_vt.0).1, fl_vt.1);
                let fl = materialize_2d_lattice(fl);

                // now collect the projected lattice
                let pl: Vec<_> = ls.iter().map(|&(x, y, _)| (((), x), y)).collect();
                let pl = LatticeCanonicalizable::canonicalize(pl);

                let pl = materialize_2d_lattice(pl);

                // now what is the rank of the intersection with t = 0?
                match fl.len() {
                    0 => {
                        // rank zero: Oscillator or glider, probably discard since we don't expect
                        // any interesting results.  Could analyze as oscillator/glider to give
                        // period and shift.

                        // TODO: care more?
                        eprintln!("   {:?}: rank zero: fl_vt = {:?}", result, fl_vt);
                    }
                    1 => {
                        // rank one: Wick of some sort.  Presumably all interesting although
                        // overpop-only connection may mean a lot of boring stuff here.  Projection
                        // of connection lattice into (x, y) plane (rather than intersection with t
                        // = 0) tells us stuff here.  If it is also one rank then we've got an
                        // oscillator wick and if it's two rank we have a (sideways) moving wick.

                        match pl.len() {
                            1 => {
                                eprintln!("   {:?}: rank one/one...", result);
                                // TODO
                            }
                            2 => {
                                eprintln!("   {:?}: rank one/two...", result);
                                // TODO
                            }
                            _ => {
                                panic!();
                            }
                        };
                    }
                    2 => {
                        // rank two: Real agar.

                        // TODO
                        eprintln!("   {:?}: rank two...", result);
                    }
                    _ => {
                        panic!();
                    }
                };
            }
        }
    }
}

fn search(lattice: (isize, isize, isize), flags: &Flags, s0: u64, results: &mut Vec<(u64, isize)>) {
    let mut prev_vec = Vec::new();
    let mut prev_map = HashMap::new();
    let mut s = s0;

    loop {
        if flags.get(s) {
            break;
        }

        if let Some(&idx) = prev_map.get(&s) {
            let &sc = (&prev_vec[idx..]).iter().min().unwrap();
            results.push((sc, (prev_vec.len() - idx) as isize));
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
                    let (x2, y2, _, _) = canon_2d(lattice, x + dx, y + dy);
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

fn compute_links(lattice: (isize, isize, isize), s0: u64) -> HashMap<(isize, isize, isize), HashSet<((isize, isize, isize), (isize, isize))>> {
    let (w, h, _sx) = lattice;
    let mut nss = HashMap::new();
    let mut links = HashMap::new();
    let mut add_link = |p1, p2, l: (isize, isize)| {
        let (lx, ly) = l;
        links.entry(p1).or_insert_with(|| HashSet::new()).insert((p2, (lx, ly)));
        links.entry(p2).or_insert_with(|| HashSet::new()).insert((p1, (-lx, -ly)));
    };
    for x in 0..w {
        for y in 0..h {
            let idx = y * w + x;
            if (s0 >> idx) & 1 == 1 {
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        let (x2, y2, lx, ly) = canon_2d(lattice, x + dx, y + dy);
                        nss.entry((x2, y2)).or_insert_with(|| Vec::new()).push(((x, y), (-lx, -ly)));
                    }
                }
            }
        }
    }
    for ((x, y), ns) in nss.into_iter() {
        let idx = y * w + x;
        let ct = ns.len();
        let living_curr = ((s0 >> idx) & 1 == 1);

        // If we're alive or we're dead and overpopulated then link up the whole neighborhood.
        if living_curr || ct >= 3 {
            for &((x1, y1), (lx1, ly1)) in &ns {
                for &((x2, y2), (lx2, ly2)) in &ns {
                    add_link((x1, y1, 0), (x2, y2, 0), (-lx1 + lx2, -ly1 + ly2));
                }
            }
        }
        // If the future is alive, then link the neighborhood to it as well.
        let living_next = match living_curr {
            true => (3 <= ct && ct <= 4),
            false => ct == 3,
        };
        if living_next {
            for &((x1, y1), (lx1, ly1)) in &ns {
                add_link((x, y, 1), (x1, y1, 0), (lx1, ly1));
            }
        }
    }

    links
}

fn canon_1d(m: isize, a: isize) -> (isize, isize) {
    let mut a = a;
    let mut la = 0;

    while a < 0 {
        a += m;
        la -= 1;
    }
    while a >= m {
        a -= m;
        la += 1;
    }

    (a, la)
}

fn canon_2d(lattice: (isize, isize, isize), x: isize, y: isize) -> (isize, isize, isize, isize) {
    let (w, h, sx) = lattice;
    let mut x = x;
    let mut y = y;
    let mut lx = 0;
    let mut ly = 0;

    while y < 0 {
        x += sx;
        y += h;
        ly -= 1;
    }
    while y >= h {
        x -= sx;
        y -= h;
        ly += 1;
    }
    while x < 0 {
        x += w;
        lx -= 1;
    }
    while x >= w {
        x -= w;
        lx += 1;
    }
    (x, y, lx, ly)
}

fn compute_lattice_links(links: &HashMap<(isize, isize, isize), HashSet<((isize, isize, isize), (isize, isize, isize))>>) -> Option<Vec<(isize, isize, isize)>> {
    let &p1 = links.keys().next().unwrap();

    let (found, r) = algo::search_lattice_links(links, p1);

    for &p in links.keys() {
        if p.2 != 0 {
            continue;
        }

        if !found.contains(&p) {
            // All our cells (in t = 0) should have been part of same connected component or we
            // discard since we should find connected components separately.
            return None;
        }
    }

    Some(r)
}

fn materialize_2d_lattice(r: (Option<(((), isize), isize)>, (Option<((), isize)>, ()))) -> Vec<(isize, isize)> {
    let (vy, r) = r;
    let (vx, _) = r;

    let mut r = Vec::new();
    if let Some((((), x), y)) = vy {
        r.push((x, y));
    }
    if let Some(((), x)) = vx {
        r.push((x, 0));
    }
    r
}
