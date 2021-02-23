#![allow(unused_parens)]

#[macro_use]
extern crate ars_ds;

use chrono::Local;
use crossbeam::queue::PopError;
use crossbeam::queue::SegQueue;
use rand::Rng;
use rand::seq::SliceRandom;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufRead;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;
use std::time::Duration;

mod flags;
mod lattice;

use flags::Flags;
use lattice::CanonicalLattice;
use lattice::LatticeCanonicalizable;

type Vec3 = (isize, isize, isize);
type Vec2 = (isize, isize);
type Vec1 = (isize,);
type Geometry2 = (Option<Vec2>, (Option<Vec1>, ()));
type Geometry3 = (Option<Vec3>, Geometry2);

fn main() {
    let mut args = std::env::args().skip(1);
    let cmd = args.next().unwrap();

    if cmd == "genl" {
        for n in args {
            let n = n.parse().unwrap();
            genl(n);
        }
        return;
    }

    if cmd == "gens" {
        let stdin = std::io::stdin();
        for line in stdin.lock().lines() {
            let line = line.unwrap();
            let parts: Vec<_> = line.split(' ').collect();
            assert_eq!(parts.len(), 3);

            let mx = parts[0].parse().unwrap();
            let my = parts[1].parse().unwrap();
            let syx = parts[2].parse().unwrap();

            gens(mx, my, syx);
        }
        return;
    }

    if cmd == "anas" {
        let stdin = std::io::stdin();
        for line in stdin.lock().lines() {
            let line = line.unwrap();
            let parts: Vec<_> = line.split(' ').collect();
            assert_eq!(parts.len(), 8);

            let mx = parts[0].parse().unwrap();
            let my = parts[1].parse().unwrap();
            let syx = parts[2].parse().unwrap();
            let s = parts[3].parse().unwrap();
            let period = parts[4].parse().unwrap();

            for result in anas(mx, my, syx, s, period) {
                print_res(&result);
            }
        }
        return;
    }

    if cmd == "rand" {
        return main_rand();
    }

    panic!("Unknown cmd {:?}", cmd);
}

fn debug_log(msg: impl AsRef<str>) {
    let msg = msg.as_ref();
    println!("{} - {}", Local::now().format("%Y%m%d %H:%M:%S"), msg);
}

fn debug_time<T>(label: impl AsRef<str>, cb: impl FnOnce() -> T) -> T {
    let label = label.as_ref();
    let t0 = std::time::Instant::now();
    debug_log(format!("Starting {}...", label));
    let ret = cb();
    debug_log(format!("Finished {}: {:?}", label, t0.elapsed()));
    return ret;
}

fn genl(n: isize) {
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
            let t_syx = (Some((t_mx,)), ()).canonicalize((t_syx,)).0;
            // if smaller flipped...
            let t_syx = t_syx.min(t_mx - t_syx);

            if (-t_mx, t_h, t_syx) < (-mx, my, syx) {
                continue;
            }

            println!("{} {} {}", mx, my, syx);
        }
    }
}

fn gens(mx: isize, my: isize, syx: isize) {
    let lattice = (mx, my, syx);
    let n = mx * my;

    let t0 = std::time::Instant::now();

    let threads = 8;
    let workunit_bits = 6.min(n);

    let geometry2 = (Some((syx, my)), (Some((mx,)), ()));

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

    eprintln!("Lattice {:?} took {:?}", lattice, t0.elapsed());
    eprintln!("Lattice {:?} => {} results", lattice, results.len());

    for (s, period, t, dx, dy) in results {
        println!("{} {} {} {} {} {} {} {}", mx, my, syx, s, period, t, dx, dy);
    }
}

fn anas(mx: isize, my: isize, syx: isize, s: u64, period: isize) -> Vec<(Geometry3, Vec<Vec2>)> {
    let mut gen0 = HashSet::new();
    for x in 0..mx {
        for y in 0..my {
            let idx = y * mx + x;
            if (s >> idx) & 1 == 1 {
                gen0.insert((x, y));
            }
        }
    }
    ana2(mx, my, syx, period, gen0)
}

fn ana2(mx: isize, my: isize, syx: isize, ttp: isize, gen0: HashSet<Vec2>) -> Vec<(Geometry3, Vec<Vec2>)> {
    let geometry3 = (Some((0, 0, ttp)), (Some((syx, my)), (Some((mx,)), ())));

    let mut links = HashMap::new();
    let mut gen1 = gen0.clone();
    for t in 0..ttp {
        let (links1, gen2) = compute_step_links(mx, my, syx, &gen1);
        for ((x1, y1, t1), links2) in links1.into_iter() {
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
        gen1 = gen2
    }
    assert!(gen0 == gen1);

    // links is the directed graph of (x, y, t) with edges labelled with absolute shift

    let mut checked = HashSet::new();
    let mut ret = Vec::new();
    for &p1 in links.keys() {
        if p1.2 != 0 {
            continue;
        }

        if checked.contains(&p1) {
            continue;
        }

        // These "lattice links" tell us the generators for the lattice that connects the same cell
        // in each tile (in the 3D lattice of space and time).  These are given in absolute
        // coordinates, not wrap counts.
        let lat = ars_graph::weighted::find_cycle_generators(&links, p1);
        let lat = LatticeCanonicalizable::canonicalize(lat);

        // The values here tell us a lattice distance by which p1 is connected to the key.  Again,
        // given in absolute coordinates.
        let connected = ars_graph::weighted::find_connected(&links, p1);

        // Mark these as handled so we don't revisit.
        for &p2 in connected.keys() {
            if p2.2 != 0 {
                continue;
            }
            checked.insert(p2);
        }

        // Collect cells in their positions in the bigger lattice (dictated by ls, not by
        // mx/my/syx)
        let cells = connected.iter().map(|(&(x2, y2, t2), &(dx, dy, dt))| lat.canonicalize((x2 + dx, y2 + dy, t2 + dt))).collect();

        ana2b(&mut ret, lat, cells);
    }

    ret
}

fn ana2b(ret: &mut Vec<(Geometry3, Vec<Vec2>)>, mut lat: Geometry3, mut cells: HashSet<Vec3>) {
    'top: loop {
        let &(x1, y1, t1) = cells.iter().next().unwrap();
        for &(x2, y2, t2) in &cells {
            if (x1, y1, t1) == (x2, y2, t2) {
                continue;
            }
            let dx = x2 - x1;
            let dy = y2 - y1;
            let dt = t2 - t1;
            let cells2 = cells.iter().map(|&(x, y, t)| lat.canonicalize((x + dx, y + dy, t + dt))).collect();
            if cells == cells2 {
                let mut lat2 = lat.materialize();
                lat2.push((dx, dy, dt));

                lat = LatticeCanonicalizable::canonicalize(lat2);
                cells = cells.into_iter().map(|p| lat.canonicalize(p)).collect();

                continue 'top;
            }
        }

        break;
    }

    let mut tups = Vec::new();
    let mut by_t = HashMap::new();
    for &(x0, y0, t0) in cells.iter() {
        by_t.entry(t0).or_insert_with(|| HashSet::new()).insert((x0, y0));
    }
    let min_pop = by_t.iter().map(|(_, cells)| cells.len()).min().unwrap();
    for cells in by_t.values() {
        if cells.len() != min_pop {
            continue;
        }
        for &(x0, y0) in cells {
            for &flip_x in &[false, true] {
                for &flip_y in &[false, true] {
                    for &swap_xy in &[false, true] {
                        let mangle = |(mut x, mut y): Vec2| {
                            if flip_x {
                                x = -x;
                            }
                            if flip_y {
                                y = -y;
                            }
                            if swap_xy {
                                core::mem::swap(&mut x, &mut y);
                            }
                            (x, y)
                        };

                        let lat1 = LatticeCanonicalizable::canonicalize(lat.materialize().into_iter().map(|(x, y, t)| {
                            let (x, y) = mangle((x, y));
                            (x, y, t)
                        }).collect());
                        let lat1_2d = lat1.1;
                        let mut cells1: Vec<Vec2> = cells.iter().map(|&(x, y)| lat1_2d.canonicalize(mangle((x - x0, y - y0)))).collect();
                        cells1.sort();
                        let (stx, sty, _) = lat1.0.unwrap();
                        tups.push((stx.abs() + sty.abs(), lat1, cells1));
                    }
                }
            }
        }
    }
    let (_, lat, cells) = tups.into_iter().min().unwrap();

    ret.push((lat, cells));
}

fn print_res(result: &(Geometry3, Vec<Vec2>)) {
    let &(lat, _) = result;
    let (stx, sty, mt) = lat.0.unwrap();
    let lat2d = lat.1;

    // now what is the rank of the intersection with t = 0?
    match lat2d.materialize().len() {
        0 => {
            // rank zero: Oscillator or glider, probably discard since we don't expect
            // any interesting results.  Could analyze as oscillator/glider to give
            // period and shift.

            if mt == 1 {
                assert_eq!(0, stx);
                assert_eq!(0, sty);
                println!("{:?}: still life", result);
            }
            else {
                if stx == 0 && sty == 0 {
                    println!("{:?}: p{} oscillator", result, mt);
                }
                else {
                    println!("{:?}: {} space ship", result, pretty_speed(stx, sty, mt));
                }
            }
        }
        1 => {
            // rank one: Wick of some sort.  Presumably all interesting although overpop-only
            // connection may mean a lot of boring stuff here.

            // Redo canonicalization with t innermost to try to figure out what's going on...
            let bw_lat = LatticeCanonicalizable::canonicalize(lat.materialize().into_iter().map(|(x, y, t)| (t, x, y)).collect());
            match ((bw_lat.1).1).0 {
                Some((ttp,)) => {
                    // We have a stationary period: wick.

                    if (stx, sty, mt) == (0, 0, 1) {
                        // Didn't really need to analyze bw_lat for this...
                        println!("{:?}: still life wick", result);
                    }
                    else if (stx, sty, mt) == (0, 0, ttp) {
                        println!("{:?}: p{} oscillator wick", result, mt);
                    }
                    else {
                        println!("{:?}: {} jump p{} oscillator wick", result, pretty_speed(stx, sty, mt), ttp);
                    }
                },
                None => {
                    // No stationary period: wave.

                    // We want to show the "smallest" possible interpretted speed ordered by
                    // movement distance, not by period.  E.g.  something whose minimum period jump
                    // is (2, 1)c might be readable as 2c/2 which is strictly better (sort for
                    // smallest total travel, break ties by period).

                    let lat_2d = lat.1;
                    // TODO: wtf is max relevant multipler?
                    let (_, mult, (stx1, sty1, mt1)) = (1..100).map(|mult| {
                        let (stx1, sty1) = lat_2d.canonicalize((mult * stx, mult * sty));
                        (stx1.abs() + sty1.abs(), mult, (stx1, sty1, mult * mt))
                    }).min().unwrap();

                    if mult == 1 {
                        println!("{:?}: {} wave", result, pretty_speed(stx, sty, mt));
                    }
                    else {
                        println!("{:?}: {} wave (jump {})", result, pretty_speed(stx1, sty1, mt1), pretty_speed(stx, sty, mt));
                    }
                },
            }
        }
        2 => {
            // rank two: Real agar.

            if (stx, sty, mt) == (0, 0, 1) {
                println!("{:?}: still life agar", result);
            }
            else if (stx, sty) == (0, 0) {
                println!("{:?}: p{} agar", result, mt);
            }
            else {
                let bw_lat = LatticeCanonicalizable::canonicalize(lat.materialize().into_iter().map(|(x, y, t)| (t, x, y)).collect());
                let (ttp,) = ((bw_lat.1).1).0.unwrap();
                println!("{:?}: jump {} p{} agar", result, pretty_speed(stx, sty, mt), ttp);
            }
        }
        _ => {
            panic!();
        }
    }
}

fn search(lattice: Vec3, masks: &Vec<Vec<u64>>, flags: &Flags, s0: u64, results: &mut HashSet<(u64, isize)>) {
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

fn tick(lattice: Vec3, masks: &Vec<Vec<u64>>, s0: u64) -> u64 {
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

fn compute_step_links(mx: isize, my: isize, syx: isize, gen0: &HashSet<Vec2>) -> (HashMap<Vec3, HashSet<(Vec3, Vec2)>>, HashSet<Vec2>) {
    let geometry2 = (Some((syx, my)), (Some((mx,)), ()));
    let mut nss = HashMap::new();
    let mut links = HashMap::new();
    let mut add_link = |p1, p2, l: Vec2| {
        let (lx, ly) = l;
        links.entry(p1).or_insert_with(|| HashSet::new()).insert((p2, (lx, ly)));
        links.entry(p2).or_insert_with(|| HashSet::new()).insert((p1, (-lx, -ly)));
    };
    for &(x, y) in gen0 {
        for dx in -1..=1 {
            for dy in -1..=1 {
                let x2 = x + dx;
                let y2 = y + dy;
                let (cx2, cy2) = geometry2.canonicalize((x2, y2));
                nss.entry((cx2, cy2)).or_insert_with(|| Vec::new()).push(((x, y), (cx2 - x2, cy2 - y2)));
            }
        }
    }
    let mut gen1 = HashSet::new();
    for ((x, y), ns) in nss.into_iter() {
        let ct = ns.len();
        let living_curr = gen0.contains(&(x, y));

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
            gen1.insert((x, y));
        }
    }

    (links, gen1)
}

fn pretty_speed(x: isize, y: isize, mt: isize) -> String {
    let x = x.abs();
    let y = y.abs();

    let (x, y) = (x.max(y), x.min(y));

    let mut ret;
    if y == 0 {
        if x == 1 {
            ret = format!("c");
        }
        else {
            ret = format!("{}c", x);
        }
    }
    else {
        ret = format!("({}, {})c", x, y);
    }

    if mt != 1 {
        ret = format!("{}/{}", ret, mt);
    }

    return ret;
}

fn compute_masks(lattice: Vec3) -> Vec<Vec<u64>> {
    let (mx, my, syx) = lattice;
    let geometry2 = (Some((syx, my)), (Some((mx,)), ()));

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

fn main_rand() {
    let threads = 8;
    let sleep_ms = 5000;
    let max_lattice = 1000;

    let mut already = HashSet::new();
    let mut lattices = HashMap::new();
    for mx in 1..=max_lattice {
        for my in 1..=(max_lattice / mx) {
            let sz = mx * my;
            let lats = lattices.entry(sz).or_insert_with(|| Vec::new());
            for syx in 0..=(mx / 2) {
                lats.push((mx, my, syx));
            }
        }
    }
    let lattices: Vec<_> = lattices.into_iter().map(|(_, lats)| lats).collect();

    loop {
        debug_time("round", || {
            let mut results: Vec<_> = (0..threads).map(|_| Vec::new()).collect();
            let stop = AtomicBool::new(false);

            {
                crossbeam::scope(|sc| {
                    for results in results.iter_mut() {
                        let already = &already;
                        let lattices = &lattices;
                        let stop = &stop;
                        sc.spawn(move |_| {
                            while !stop.load(Ordering::Relaxed) {
                                for result in main_rand1(stop, lattices) {
                                    if already.contains(&result) {
                                        continue;
                                    }
                                    results.push(result);
                                }
                            }
                        });
                    }

                    std::thread::sleep(Duration::from_millis(sleep_ms));

                    stop.store(true, Ordering::Relaxed);
                }).unwrap();
            }

            for results in results.into_iter() {
                for result in results.into_iter() {
                    if already.contains(&result) {
                        continue
                    }
                    print_res(&result);
                    already.insert(result);
                }
            }
        });
    }
}

fn compute_step(mx: isize, my: isize, syx: isize, gen0: &HashSet<Vec2>) -> HashSet<Vec2> {
    let geometry2 = (Some((syx, my)), (Some((mx,)), ()));
    let mut cts = HashMap::new();
    for &(x, y) in gen0 {
        for dx in -1..=1 {
            for dy in -1..=1 {
                let x2 = x + dx;
                let y2 = y + dy;
                let (cx2, cy2) = geometry2.canonicalize((x2, y2));
                *cts.entry((cx2, cy2)).or_insert_with(|| 0) += 1;
            }
        }
    }
    let mut gen1 = HashSet::new();
    for ((x, y), ct) in cts.into_iter() {
        let living_curr = gen0.contains(&(x, y));

        // If the future is alive, then link the neighborhood to it as well.
        let living_next = match living_curr {
            true => (3 <= ct && ct <= 4),
            false => ct == 3,
        };
        if living_next {
            gen1.insert((x, y));
        }
    }

    gen1
}

fn main_rand1(stop: &AtomicBool, lattices: &[Vec<Vec3>]) -> Vec<(Geometry3, Vec<Vec2>)> {
    let lats = lattices.choose(&mut rand::thread_rng()).unwrap();
    let &(mx, my, syx) = lats.choose(&mut rand::thread_rng()).unwrap();

    let mut gen0 = Vec::new();
    for x in 0..mx {
        for y in 0..my {
            if rand::thread_rng().gen() {
                gen0.push((x, y));
            }
        }
    }
    gen0.sort();

    let mut t = 0isize;
    let mut already = HashMap::new();
    loop {
        if stop.load(Ordering::Relaxed) {
            return vec![];
        }
        if let Some(t0) = already.get(&gen0) {
            if t - t0 > 200 {
                return vec![];
            }
            return ana2(mx, my, syx, t - t0, gen0.iter().cloned().collect());
        }
        already.insert(gen0.clone(), t);

        let gen1 = compute_step(mx, my, syx, &gen0.iter().cloned().collect());
        let mut gen1: Vec<_> = gen1.into_iter().collect();
        gen1.sort();

        t += 1;
        gen0 = gen1;
    }
}
