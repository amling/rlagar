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
mod tuple;

use flags::Flags;
use lattice::Canonicalizes;
use lattice::LatticeCanonicalizable;

fn main() {
for n in 2..=23 {
let t0 = std::time::Instant::now();
    let threads = 8;
    let workunit_bits = 6.min(n);

    let mut lattices = Vec::new();
    for mx in 1..=n {
        if n % mx != 0 {
            continue;
        }

        let my = n / mx;

        // syx is the shift in the x direction when we wrap in the y direction
        for syx in 0..mx {
            if 2 * syx > mx {
                // don't bother with shifts over half, instead we'll find whatever flipped
                continue;
            }

            let (d, _, (s, t)) = lattice::egcd(syx, mx, (1, 0), (0, 1));
            // d = s * syx + t * mx

            // (mx, 0), (syx, my)
            // (mx, 0), (syx, my), (syx * mx / d, my * mx / d), (s * syx + t * mx, s * my + t * 0)
            // (mx, 0), (syx, my), (syx/d * mx, my * mx / d), (d, s * my)
            // (mx, 0), (syx, my), (0, my * mx / d), (d, s * my)

            // the order stepping vertically, what mx would be if we transposed
            let t_mx = my * mx / d;
            // the rest...
            let ht = d;
            // what syx would be if we transposed
            let t_syx = s * t;
            // grumble grumble stupid mod
            let t_syx = lattice::unvec1(lattice::geom1(t_mx).canonicalize(lattice::vec1(t_syx)));
            // if smaller flipped...
            let t_syx = t_syx.min(t_mx - t_syx);

            if (-t_mx, ht, t_syx) < (-mx, ht, syx) {
                continue;
            }

            lattices.push((mx, my, syx));
        }
    }
    //println!("Lattices {:?}", lattices);

    for lattice in lattices {
        let (mx, my, syx) = lattice;
        let geometry2 = lattice::geom2(mx, my, syx);

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

        let results = results.filter_map(|(s, period)| {
            let mut s1 = s;
            for t in 0.. {
                for dx in 0..mx {
                    for dy in 0..my {
                        if dx == 0 && dy == 0 && t == 0 {
                            continue;
                        }

                        let mut s2 = 0;
                        for x in 0..mx {
                            for y in 0..my {
                                let idx = y * mx + x;

                                let (x2, y2) = lattice::unvec2(geometry2.canonicalize(lattice::vec2(x + dx, y + dy)));
                                let idx2 = y2 * mx + x2;

                                if (s1 >> idx) & 1 == 1 {
                                    s2 |= (1 << idx2);
                                }
                            }
                        }

                        if s2 < s {
                            // found lower value somewhere else, drop
                            //eprintln!("Drop lattice {:?} result {}, shift ({}, {}) at t {}", lattice, s, dx, dy, t);
                            return None;
                        }

                        if s2 == s {
                            if t == 0 {
                                // found extra spatial symmetry, this can be found in smaller
                                // geometry
                                //eprintln!("Drop lattice {:?} result {} shift ({}, {})", lattice, s, dx, dy);
                                return None;
                            }
                            let (dx, dy) = lattice::unvec2(geometry2.canonicalize(lattice::vec2(-dx, -dy)));
                            return Some((s, period, t, dx, dy));
                        }
                    }
                }

                s1 = tick(lattice, s1);
            }
            panic!();
        });

        let results: BTreeSet<_> = results.collect();
        eprintln!("Lattice {:?} => {} results", lattice, results.len());
        for result in results {
            //eprintln!("   {:?}", result);

            let (s, period, mt, stx, sty) = result;
            let geometry3 = lattice::geom3(mx, my, syx, mt, stx, sty);

            let mut links = HashMap::new();
            let mut s1 = s;
            for t in 0..mt {
                for ((x1, y1, t1), links2) in compute_links(lattice, s1).into_iter() {
                    let t1 = t + t1;
                    let (cx1, cy1, ct1) = lattice::unvec3(geometry3.canonicalize(lattice::vec3(x1, y1, t1)));
                    for ((x2, y2, t2), (dx, dy)) in links2.into_iter() {
                        let t2 = t + t2;
                        let (cx2, cy2, ct2) = lattice::unvec3(geometry3.canonicalize(lattice::vec3(x2, y2, t2)));
                        links.entry((cx1, cy1, ct1)).or_insert_with(|| HashSet::new()).insert(((cx2, cy2, ct2), (x2 - x1 + dx, y2 - y1 + dy, t2 - t1)));
                    }
                }
                s1 = tick(lattice, s1);
            }

            // links is the directed graph of (x, y, t) with edges labelled with absolute shift

            if let Some(ls) = compute_lattice_links(&links) {
                // These "lattice links" tell us the generators for the lattice that connects the
                // same cell in each fundamntal period (in the 3D lattice of space and time).
                // These are given in absolute x/y coordinates, not wrap counts.

                let fl: Vec<_> = ls.iter().map(|&(x, y, t)| lattice::vec3(x, y, t)).collect();
                let fl = LatticeCanonicalizable::canonicalize(fl);

                let (fl_vt, fl) = fl;
                let fl_vt = lattice::unvec3(fl_vt.unwrap());
                let fl = materialize_2d_lattice(fl);

                // now collect the projected lattice
                let pl: Vec<_> = ls.iter().map(|&(x, y, _)| lattice::vec2(x, y)).collect();
                let pl = LatticeCanonicalizable::canonicalize(pl);

                let pl = materialize_2d_lattice(pl);

                // now what is the rank of the intersection with t = 0?
                match fl.len() {
                    0 => {
                        // rank zero: Oscillator or glider, probably discard since we don't expect
                        // any interesting results.  Could analyze as oscillator/glider to give
                        // period and shift.

                        let (stx, sty, mt) = fl_vt;
                        if mt == 1 {
                            assert_eq!(0, stx);
                            assert_eq!(0, sty);
                            eprintln!("   {}: still life", s);
                        }
                        else {
                            if stx == 0 && sty == 0 {
                                assert_eq!(period, mt);
                                eprintln!("   {}: oscillator, period {}", s, mt);
                            }
                            else {
                                eprintln!("   {}: space ship, shift ({}, {}), period {}", s, stx, sty, mt);
                            }
                        }
                    }
                    1 => {
                        // rank one: Wick of some sort.  Presumably all interesting although
                        // overpop-only connection may mean a lot of boring stuff here.  Projection
                        // of connection lattice into (x, y) plane (rather than intersection with t
                        // = 0) tells us stuff here.  If it is also one rank then we've got an
                        // oscillator wick and if it's two rank we have a (sideways) moving wick.

                        match pl.len() {
                            1 => {
                                let (stx, sty, mt) = fl_vt;
                                if stx == 0 && sty == 0 && mt == 1 {
                                    eprintln!("   {}: still life wick", s);
                                }
                                else {
                                    if stx == 0 && sty == 0 {
                                        assert_eq!(period, mt);
                                        eprintln!("   {}: oscillator wick, period {}", s, mt);
                                    }
                                    else {
                                        // not actual period when including a shift
                                        eprintln!("   {}: shifting oscillator wick, shift ({}, {}), period {}, true period {}", s, stx, sty, mt, period);
                                    }
                                }
                            }
                            2 => {
                                // TODO
                                eprintln!("   {:?}: space ship wick", result);
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
eprintln!("N {} took {:?}", n, t0.elapsed());
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
    let (mx, my, sx) = lattice;
    let geometry2 = lattice::geom2(mx, my, sx);
    let mut s1 = 0;
    for x in 0..mx {
        for y in 0..my {
            let idx = y * mx + x;
            let mut s = 0;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    let (x2, y2) = lattice::unvec2(geometry2.canonicalize(lattice::vec2(x + dx, y + dy)));
                    let idx2 = y2 * mx + x2;
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
    let (mx, my, sx) = lattice;
    let geometry2 = lattice::geom2(mx, my, sx);
    let mut nss = HashMap::new();
    let mut links = HashMap::new();
    let mut add_link = |p1, p2, l: (isize, isize)| {
        let (lx, ly) = l;
        links.entry(p1).or_insert_with(|| HashSet::new()).insert((p2, (lx, ly)));
        links.entry(p2).or_insert_with(|| HashSet::new()).insert((p1, (-lx, -ly)));
    };
    for x in 0..mx {
        for y in 0..my {
            let idx = y * mx + x;
            if (s0 >> idx) & 1 == 1 {
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        let x2 = x + dx;
                        let y2 = y + dy;
                        let (cx2, cy2) = lattice::unvec2(geometry2.canonicalize(lattice::vec2(x2, y2)));
                        nss.entry((cx2, cy2)).or_insert_with(|| Vec::new()).push(((x, y), (cx2 - x2, cy2 - y2)));
                    }
                }
            }
        }
    }
    for ((x, y), ns) in nss.into_iter() {
        let idx = y * mx + x;
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
    if let Some(vy) = vy {
        r.push(lattice::unvec2(vy));
    }
    if let Some(vx) = vx {
        r.push((lattice::unvec1(vx), 0));
    }
    r
}
