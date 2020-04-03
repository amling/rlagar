#![allow(unused_parens)]

#[macro_use]
extern crate ars_ds;

use ars_ds::tuple::Tuple1;
use crossbeam::queue::PopError;
use crossbeam::queue::SegQueue;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;

mod flags;
mod lattice;

use flags::Flags;
use lattice::Canonicalizes;
use lattice::LatticeCanonicalizable;

fn main() {
for n in 2..=33 {
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

            let (d, _, (s, _)) = ars_aa::misc::egcd(syx, mx, (1, 0), (0, 1));
            // d = s * syx + t * mx

            // (mx, 0), (syx, my)
            // (mx, 0), (syx, my), (syx * mx / d, my * mx / d), (s * syx + t * mx, s * my + t * 0)
            // (mx, 0), (syx, my), (syx/d * mx, my * mx / d), (d, s * my)
            // (mx, 0), (syx, my), (0, my * mx / d), (d, s * my)

            // the order stepping vertically, what mx would be if we transposed
            let t_mx = my * mx / d;
            // the rest...
            let t_h = d;
            // what syx would be if we transposed
            let t_syx = s * my;
            // grumble grumble stupid mod
            let t_syx = (Some(Tuple1(t_mx)), ()).canonicalize(Tuple1(t_syx)).0;
            // if smaller flipped...
            let t_syx = t_syx.min(t_mx - t_syx);

            if (-t_mx, t_h, t_syx) < (-mx, my, syx) {
                continue;
            }

            lattices.push((mx, my, syx));
        }
    }
    //println!("Lattices {:?}", lattices);

    for lattice in lattices {
        let (mx, my, syx) = lattice;
        let geometry2 = (Some((syx, my)), (Some(Tuple1(mx)), ()));

        let masks = compute_masks(lattice);
        let masks = &masks;

        let flags = Flags::new(1 << n);
        let flags = &flags;

        let workunits: Vec<_> = (0..(1 << workunit_bits)).collect();
        let mut results: Vec<_> = workunits.iter().map(|_| HashSet::new()).collect();

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
                                search(lattice, masks, flags, s0, results);
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

                                let (x2, y2) = geometry2.canonicalize((x + dx, y + dy));
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
                            let (dx, dy) = geometry2.canonicalize((-dx, -dy));
                            return Some((s, period, t, dx, dy));
                        }
                    }
                }

                s1 = tick(lattice, masks, s1);
            }
            panic!();
        });

        let results: BTreeSet<_> = results.collect();
        eprintln!("Lattice {:?} => {} results", lattice, results.len());
        'result: for result in results {
            //eprintln!("   {:?}", result);

            let (s, period, mt, stx, sty) = result;
            let geometry3 = (Some((stx, sty, mt)), (Some((syx, my)), (Some(Tuple1(mx)), ())));

            let mut links = HashMap::new();
            let mut s1 = s;
            for t in 0..mt {
                for ((x1, y1, t1), links2) in compute_links(lattice, s1).into_iter() {
                    let t1 = t + t1;
                    let (cx1, cy1, ct1) = geometry3.canonicalize((x1, y1, t1));
                    for ((x2, y2, t2), (lx, ly)) in links2.into_iter() {
                        let t2 = t + t2;
                        let (cx2, cy2, ct2) = geometry3.canonicalize((x2, y2, t2));

                        let rx = (cx1 - x1) + lx + (x2 - cx2);
                        let ry = (cy1 - y1) + ly + (y2 - cy2);
                        let rt = (ct1 - t1) + (t2 - ct2);

                        links.entry((cx1, cy1, ct1)).or_insert_with(|| HashSet::new()).insert(((cx2, cy2, ct2), (rx, ry, rt)));
                    }
                }
                s1 = tick(lattice, masks, s1);
            }

            // links is the directed graph of (x, y, t) with edges labelled with absolute shift

            // min not necessary for algo in general but makes debugging less insane
            let &p1 = links.keys().filter(|p| p.2 == 0).min().unwrap();

            // These "lattice links" tell us the generators for the lattice that connects the same
            // cell in each fundamental period (in the 3D lattice of space and time).  These are
            // given in absolute x/y coordinates, not wrap counts.
            let ls = ars_graph::weighted::find_cycle_generators(&links, p1);
            let fl = ls.clone();
            let fl = LatticeCanonicalizable::canonicalize(fl);

            // The values here tell us a lattice distance by which p1 is connected to the key.
            // Again, given in absolute coordinates.
            let connected = ars_graph::weighted::find_connected(&links, p1);

            for &p2 in links.keys() {
                if p2.2 != 0 {
                    continue;
                }

                // All our cells (in t = 0) should have been part of same connected component or we
                // discard since we should find connected components separately.
                match connected.get(&p2) {
                    Some(&pd) => {
                        let pd = fl.canonicalize(pd);
                        // It's okay if we don't connect at (0, 0, 0), we just need to connect at
                        // some point in space at time 0.
                        if pd.2 != 0 {
                            // This is pretty amazing.  This happens when e.g.  you have two copies
                            // of the same wick that appear to replace each other.  Each cell in
                            // space is connected, but not all back at t = 0 so things actually can
                            // be split into pieces.
                            continue 'result;
                        }
                    }
                    None => {
                        // not connected to this cell at all!
                        continue 'result;
                    }
                }
            }

            let (fl_vt, fl2) = fl;
            let fl_vt = fl_vt.unwrap();

            // now collect the projected lattice
            let pl: Vec<_> = ls.iter().map(|&(x, y, _)| (x, y)).collect();
            let pl = LatticeCanonicalizable::canonicalize(pl);

            // now what is the rank of the intersection with t = 0?
            match materialize_2d_lattice(fl2).len() {
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
                            eprintln!("   {}: p{} oscillator", s, mt);
                        }
                        else {
                            eprintln!("   {}: {} space ship", s, pretty_speed(fl2, mt, stx, sty));
                        }
                    }
                }
                1 => {
                    // rank one: Wick of some sort.  Presumably all interesting although
                    // overpop-only connection may mean a lot of boring stuff here.  Projection
                    // of connection lattice into (x, y) plane (rather than intersection with t
                    // = 0) tells us stuff here.  If it is also one rank then we've got an
                    // oscillator wick and if it's two rank we have a (sideways) moving wick.

                    match materialize_2d_lattice(pl).len() {
                        1 => {
                            let (stx, sty, mt) = fl_vt;
                            if stx == 0 && sty == 0 && mt == 1 {
                                eprintln!("   {}: still life wick", s);
                            }
                            else {
                                if stx == 0 && sty == 0 {
                                    assert_eq!(period, mt);
                                    eprintln!("   {}: p{} oscillator wick", s, mt);
                                }
                                else {
                                    // not actual period when including a shift
                                    eprintln!("   {}: {} shifting oscillator wick (true period {})", s, pretty_speed(fl2, mt, stx, sty), period);
                                }
                            }
                        }
                        2 => {
                            // TODO
                            let (stx, sty, mt) = fl_vt;
                            let canonical = ssw_canonical(&links, fl);
                            eprintln!("   {:?}: {} space ship wick, canonical {:?}", result, pretty_speed(fl2, mt, stx, sty), canonical);
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
            }
        }
    }
eprintln!("N {} took {:?}", n, t0.elapsed());
}
}

fn search(lattice: (isize, isize, isize), masks: &Vec<Vec<u64>>, flags: &Flags, s0: u64, results: &mut HashSet<(u64, isize)>) {
    let mut prev_vec = Vec::new();
    let mut prev_map = HashMap::new();
    let mut s = s0;

    loop {
        if flags.get(s) {
            break;
        }

        if let Some(&idx) = prev_map.get(&s) {
            let &sc = (&prev_vec[idx..]).iter().min().unwrap();
            results.insert((sc, (prev_vec.len() - idx) as isize));
            break;
        }

        let t = prev_vec.len();
        prev_vec.push(s);
        prev_map.insert(s, t);
        s = tick(lattice, masks, s);
    }

    for sp in prev_vec {
        flags.set(sp);
    }
}

fn tick(lattice: (isize, isize, isize), masks: &Vec<Vec<u64>>, s0: u64) -> u64 {
    let (mx, my, _syx) = lattice;
    let mut s1 = 0;
    for idx in 0..(mx * my) {
        let mut s = 0;
        for mask in masks[idx as usize].iter() {
            s += (s0 & mask).count_ones();
        }
        let living = match ((s0 >> idx) & 1 == 1) {
            true => (2 <= s && s <= 3),
            false => (s == 3),
        };
        if living {
            s1 |= (1 << idx);
        }
    }
    s1
}

fn compute_links(lattice: (isize, isize, isize), s0: u64) -> HashMap<(isize, isize, isize), HashSet<((isize, isize, isize), (isize, isize))>> {
    let (mx, my, syx) = lattice;
    let geometry2 = (Some((syx, my)), (Some(Tuple1(mx)), ()));
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
                        let (cx2, cy2) = geometry2.canonicalize((x2, y2));
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

fn materialize_2d_lattice(r: (Option<(isize, isize)>, (Option<Tuple1<isize>>, ()))) -> Vec<(isize, isize)> {
    let (vy, r) = r;
    let (vx, _) = r;

    let mut r = Vec::new();
    if let Some(vy) = vy {
        r.push(vy);
    }
    if let Some(Tuple1(vx)) = vx {
        r.push((vx, 0));
    }
    r
}

fn pretty_speed(geometry2: (Option<(isize, isize)>, (Option<Tuple1<isize>>, ())), mt: isize, stx: isize, sty: isize) -> String {
    // Sigh, not obvious how to make this less stupid, but it should be fine for how little it's
    // used.
    for d in 0..100 {
        for x in 0..=d {
            let y = d - x;
            for &sx in &[-1, 1] {
                for &sy in &[-1, 1] {
                    if geometry2.canonicalize((sx * x, sy * y)) == (stx, sty) {
                        let x = x.abs();
                        let y = y.abs();

                        let (x, y) = (x.max(y), x.min(y));

                        if x == 0 {
                            return "0".to_string();
                        }
                        if y == 0 {
                            if x == 1 {
                                return format!("c/{}", mt);
                            }
                            return format!("{}c/{}", x, mt);
                        }
                        return format!("({}, {})c/{}", x, y, mt);
                    }
                }
            }
        }
    }

    panic!();
}

fn compute_masks(lattice: (isize, isize, isize)) -> Vec<Vec<u64>> {
    let (mx, my, syx) = lattice;
    let geometry2 = (Some((syx, my)), (Some(Tuple1(mx)), ()));

    let mut acc = Vec::new();
    for idx in 0..(mx * my) {
        let x = idx % mx;
        let y = idx / mx;
        let mut masks = Vec::new();
        let mut add_mask = |idx2| {
            let m = 1u64 << idx2;
            for mask in masks.iter_mut() {
                if *mask & m == 0 {
                    *mask |= m;
                    return;
                }
            }
            masks.push(0);
            *masks.last_mut().unwrap() |= m;
        };
        for dx in -1..=1 {
            for dy in -1..=1 {
                if (dx, dy) == (0, 0) {
                    continue;
                }
                let (x2, y2) = geometry2.canonicalize((x + dx, y + dy));
                let idx2 = y2 * mx + x2;
                add_mask(idx2);
            }
        }
        acc.push(masks);
    }
    acc
}

// This still isn't right because once we've fixed fl flipping and swapping will not necessarily
// produce what it would have had the lattice been canonicalized that way.  Hopefully few enough
// duplicates to not go crazy...
fn ssw_canonical(links: &HashMap<(isize, isize, isize), HashSet<((isize, isize, isize), (isize, isize, isize))>>, fl: (Option<(isize, isize, isize)>, (Option<(isize, isize)>, (Option<Tuple1<isize>>, ())))) -> (usize, isize, isize, String) {
    let mut candidates = Vec::new();
    for &p1 in links.keys() {
        let connected = ars_graph::weighted::find_connected(&links, p1);
        let cells: HashSet<_> = connected.into_iter().filter_map(|(p2, pd)| {
            let pd = fl.canonicalize(pd);
            if p1.2 != p2.2 {
                return None;
            }
            assert_eq!(0, pd.2);
            Some((p2.0 - p1.0 + pd.0, p2.1 - p1.1 + pd.1))
        }).collect();
        for &mx in &[-1, 1] {
            for &my in &[-1, 1] {
                for &swap in &[false, true] {
                    let cells: HashSet<_> = cells.iter().map(|&(x, y)| {
                        let x = x * mx;
                        let y = y * my;
                        match swap {
                            true => (y, x),
                            false => (x, y),
                        }
                    }).collect();

                    let min_x = cells.iter().map(|&(x, _)| x).min().unwrap();
                    let max_x = cells.iter().map(|&(x, _)| x).max().unwrap();
                    let min_y = cells.iter().map(|&(_, y)| y).min().unwrap();
                    let max_y = cells.iter().map(|&(_, y)| y).max().unwrap();

                    let mut cells_str = String::new();
                    for y in min_y..=max_y {
                        for x in min_x..=max_x {
                            cells_str.push(match cells.contains(&(x, y)) {
                                true => '*',
                                false => '.',
                            });
                        }
                        if y != max_y {
                            cells_str.push('/');
                        }
                    }

                    candidates.push((cells.len(), max_y - min_y + 1, max_x - min_x + 1, cells_str));
                }
            }
        }
    }
    candidates.into_iter().min().unwrap()
}
