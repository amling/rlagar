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
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::sync::Arc;
use std::sync::Condvar;
use std::sync::Mutex;
use std::time::Duration;
use std::time::Instant;

mod flags;
mod lattice;

use flags::Flags;
use flags::HackFlags;
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
        let max_mask_bits = 20;
        let fsz = args.next().unwrap().parse().unwrap();

        let stdin = std::io::stdin();
        for line in stdin.lock().lines() {
            let line = line.unwrap();
            let parts: Vec<_> = line.split(' ').collect();
            assert_eq!(parts.len(), 3);

            let mx = parts[0].parse().unwrap();
            let my = parts[1].parse().unwrap();
            let syx = parts[2].parse().unwrap();

            gens(mx, my, syx, max_mask_bits, fsz);
        }
        return;
    }

    if cmd == "rand" {
        let min_area = args.next().unwrap().parse().unwrap();
        let max_area = args.next().unwrap().parse().unwrap();
        return main_rand(min_area, max_area);
    }

    if cmd == "printr" {
        let mut already = HashSet::new();

        let mut handle_line = |line: Result<String, std::io::Error>| {
            let line = line.unwrap();
            let res = serde_json::from_str(&line).unwrap();
            if already.contains(&res) {
                return;
            }
            print_res(&res);
            already.insert(res);
        };

        let mut handle_arg = |arg: &str| {
            if arg == "-" {
                let stdin = std::io::stdin();
                for line in stdin.lock().lines() {
                    handle_line(line);
                }
            }
            else {
                let f = File::open(arg).unwrap();
                let r = BufReader::new(f);
                for line in r.lines() {
                    handle_line(line);
                }
            }
        };

        let args: Vec<_> = args.collect();
        if args.is_empty() {
            handle_arg("-");
        }
        else {
            for arg in args {
                println!("Starting {:?}...", arg);
                handle_arg(&arg);
            }
        }
        return;
    }

    panic!("Unknown cmd {:?}", cmd);
}

fn debug_log(msg: impl AsRef<str>) {
    let msg = msg.as_ref();
    eprintln!("{} - {}", Local::now().format("%Y%m%d %H:%M:%S"), msg);
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

fn gens(mx: isize, my: isize, syx: isize, max_mask_bits: usize, fsz: usize) {
    let masks = debug_time("compute_masks", || compute_masks((mx, my, syx), max_mask_bits));

    let lattice = (mx, my, syx);
    debug_time(format!("lattice {:?}", lattice), || {
        let n = mx * my;

        let threads = 8;
        let workunit_bits = 20.min(n);

        let geometry2 = (Some((syx, my)), (Some((mx,)), ()));

        let flags = HackFlags::new(fsz);
        let flags = &flags;

        let workunits: Vec<_> = (0..(1 << workunit_bits)).collect();
        let mut results: Vec<_> = workunits.iter().map(|_| HashSet::new()).collect();

        {
            let cdl = Arc::new((Condvar::new(), Mutex::new(threads)));

            let q = SegQueue::new();
            for tuple in workunits.into_iter().zip(results.iter_mut()) {
                q.push(tuple);
            }
            let total = q.len();

            crossbeam::scope(|sc| {
                for _ in 0..threads {
                    sc.spawn(|_| {
                        loop {
                            let (workunit, results) = match q.pop() {
                                Result::Ok(tuple) => tuple,
                                Result::Err(PopError) => {
                                    break;
                                }
                            };

                            let suffix_bits = n - workunit_bits;
                            for suffix in 0..(1 << suffix_bits) {
                                let s0 = (workunit << suffix_bits) | suffix;
                                search(&masks, flags, s0, results);
                            }
                        }

                        {
                            let (c, m) = &*cdl;
                            let mut m = m.lock().unwrap();
                            let n = *m - 1;
                            *m = n;
                            if n == 0 {
                                c.notify_one();
                            }
                        }
                    });
                }

                let mut last: Option<(Instant, usize)> = None;
                let (c, m) = &*cdl;
                let mut m = m.lock().unwrap();
                while *m > 0 {
                    m = c.wait_timeout(m, Duration::from_millis(60000)).unwrap().0;

                    let now = (Instant::now(), q.len());

                    if let Some(prev) = last {
                        if prev.1 != now.1 {
                            let left = (now.0 - prev.0).mul_f64(now.1 as f64).div_f64((prev.1 - now.1) as f64);
                            debug_log(format!("Completed {}/{}, estimated {:?} left...", (total - now.1), total, left));
                        }
                    }

                    last = Some(now);
                }
            }).unwrap();
        }

        let results = results.into_iter();

        // combine results from each thread and dedupe
        let results = results.flatten().collect::<HashSet<_>>().into_iter();

        // Eliminate more symmetry crap.  In particular things which have a lower-valued shift of
        // themselves.  Also things which have a tighter spatial period (which we'll have found in
        // a smaller lattice anyway).
        let results = results.filter(|&(s, period)| {
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
                            return false;
                        }

                        if s2 == s {
                            if t == 0 {
                                // found extra spatial symmetry, this can be found in smaller
                                // geometry
                                //eprintln!("Drop lattice {:?} result {} shift ({}, {})", lattice, s, dx, dy);
                                return false;
                            }
                            assert!(period % t == 0);
                            return true;
                        }
                    }
                }

                s1 = tick(&masks, s1);
            }
            panic!();
        });

        // Now actually do expensive analysis of splitting into pieces, etc.
        let results = results.flat_map(|(s, period)| {
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
        });

        let results: BTreeSet<_> = results.collect();

        debug_log(format!("Lattice {:?} => {} results", lattice, results.len()));

        for res in results {
            println!("{}", serde_json::to_string(&res).unwrap());
        }
    });
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
        'p2: for &(x2, y2, t2) in &cells {
            if (x1, y1, t1) == (x2, y2, t2) {
                continue;
            }
            let dx = x2 - x1;
            let dy = y2 - y1;
            let dt = t2 - t1;

            for &(x, y, t) in cells.iter() {
                if !cells.contains(&lat.canonicalize((x + dx, y + dy, t + dt))) {
                    continue 'p2;
                }
            }

            let mut lat2 = lat.materialize();
            lat2.push((dx, dy, dt));

            lat = LatticeCanonicalizable::canonicalize(lat2);
            cells = cells.into_iter().map(|p| lat.canonicalize(p)).collect();

            continue 'top;
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

    let print = |s| {
        println!("{:?}: {}", result, s);
    };

    let mut shifts = Vec::new();
    let show_cells;

    // now what is the rank of the intersection with t = 0?
    match &lat2d.materialize()[..] {
        [] => {
            // rank zero: Oscillator or glider, probably discard since we don't expect
            // any interesting results.  Could analyze as oscillator/glider to give
            // period and shift.

            if mt == 1 {
                assert_eq!(0, stx);
                assert_eq!(0, sty);
                print(format!("still life"));
                show_cells = false;
            }
            else {
                if stx == 0 && sty == 0 {
                    print(format!("p{} oscillator", mt));
                    show_cells = (mt > 2);
                }
                else {
                    print(format!("{} space ship", pretty_speed(stx, sty, mt)));
                    show_cells = true;
                }
            }

            shifts.push((0, 0));
        }
        [(x1, y1)] => {
            // rank one: Wick of some sort.  Presumably all interesting although overpop-only
            // connection may mean a lot of boring stuff here.

            // Redo canonicalization with t innermost to try to figure out what's going on...
            let bw_lat = LatticeCanonicalizable::canonicalize(lat.materialize().into_iter().map(|(x, y, t)| (t, x, y)).collect());
            match ((bw_lat.1).1).0 {
                Some((ttp,)) => {
                    // We have a stationary period: wick.

                    if (stx, sty, mt) == (0, 0, 1) {
                        // Didn't really need to analyze bw_lat for this...
                        print(format!("still life wick"));
                        show_cells = false;
                    }
                    else if (stx, sty, mt) == (0, 0, ttp) {
                        print(format!("p{} oscillator wick", mt));
                        show_cells = (mt > 2);
                    }
                    else {
                        print(format!("{} jump p{} oscillator wick", pretty_speed(stx, sty, mt), ttp));
                        show_cells = true;
                    }
                },
                None => {
                    // No stationary period: wave.
                    print(format!("{} wave", pretty_speed(stx, sty, mt)));
                    show_cells = true;
                },
            }

            for n in -5..=5 {
                shifts.push((n * x1, n * y1));
            }
        }
        [(x1, y1), (x2, y2)] => {
            // rank two: Real agar.

            if (stx, sty, mt) == (0, 0, 1) {
                print(format!("still life agar"));
                show_cells = false;
            }
            else if (stx, sty) == (0, 0) {
                print(format!("p{} agar", mt));
                show_cells = false;
            }
            else {
                let bw_lat = LatticeCanonicalizable::canonicalize(lat.materialize().into_iter().map(|(x, y, t)| (t, x, y)).collect());
                let (ttp,) = ((bw_lat.1).1).0.unwrap();
                print(format!("jump {} p{} agar", pretty_speed(stx, sty, mt), ttp));
                show_cells = false;
            }

            for n in -5..=5 {
                for m in -5..=5 {
                    shifts.push((n * x1 + m * x2, n * y1 + m * y2));
                }
            }
        }
        _ => {
            panic!();
        }
    }

    if show_cells {
        let mut cells = HashSet::new();
        for (dx, dy) in shifts {
            for &(x, y) in result.1.iter() {
                cells.insert((x + dx, y + dy));
            }
        }

        let min_x = cells.iter().map(|&(x, _y)| x).min().unwrap();
        let max_x = cells.iter().map(|&(x, _y)| x).max().unwrap();
        let min_y = cells.iter().map(|&(_x, y)| y).min().unwrap();
        let max_y = cells.iter().map(|&(_x, y)| y).max().unwrap();

        for y in min_y..=max_y {
            let mut r = "   ".to_string();
            for x in min_x..=max_x {
                r.push(match cells.contains(&(x, y)) {
                    true => '*',
                    false => '.',
                });
            }
            println!("{}", r)
        }
    }
}

fn search(masks: &Vec<(u64, Vec<u64>)>, flags: &impl Flags, s0: u64, results: &mut HashSet<(u64, isize)>) {
    // We unroll the first few iterations for heavy [mis]optimization.  In particular we don't
    // detect low [fundamental] period in some cases, but misdetecting them as a multiple will be
    // cleaned up later.

    if flags.get(s0) {
        return;
    }

    let s1 = tick(masks, s0);
    if flags.get(s1) {
        flags.set(s0);
        return;
    }

    let s2 = tick(masks, s1);
    if flags.get(s2) {
        flags.set(s1);
        flags.set(s0);
        return;
    }

    let mut prev_vec = vec![s0, s1, s2];
    let mut prev_map: HashMap<_, _> = vec![(s0, 0), (s1, 1), (s2, 2)].into_iter().collect();
    let mut s = s2;

    loop {
        s = tick(masks, s);

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
    }

    for sp in prev_vec {
        flags.set(sp);
    }
}

fn tick(masks: &Vec<(u64, Vec<u64>)>, s0: u64) -> u64 {
    let mut s1 = 0;
    for &(mask, ref results) in masks {
        let masked = bitintr::Pext::pext(s0, mask);
        s1 |= results[masked as usize];
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

fn compute_masks(lattice: Vec3, max_mask_bits: usize) -> Vec<(u64, Vec<u64>)> {
    let (mx, my, syx) = lattice;
    let geometry2 = (Some((syx, my)), (Some((mx,)), ()));

    let neighbors_of = |(x, y)| {
        (-1..=1).flat_map(move |dx| (-1..=1).map(move |dy| geometry2.canonicalize((x + dx, y + dy))))
    };

    let compute_masks1 = |mask_ps: HashSet<Vec2>, covered: Vec<Vec2>| {
        // Put the bits in order
        let mask_ps = {
            let mut v: Vec<_> = mask_ps.into_iter().collect();
            v.sort_by_key(|&(x, y)| y * mx + x);
            v
        };

        // Assemble our mask
        let mask = {
            let mut r = 0u64;
            for &(x, y) in mask_ps.iter() {
                let idx = y * mx + x;
                r |= (1 << idx);
            }
            r
        };

        let mut map = vec![0; 1 << mask_ps.len()];
        for masked in 0..(1 << mask_ps.len()) {
            let mut p_live = HashMap::new();
            for (i, &(x, y)) in mask_ps.iter().enumerate() {
                let live = (masked & (1 << i)) != 0;
                p_live.insert((x, y), live);
            }

            let mut res = 0;
            for &(x, y) in covered.iter() {
                let mut ct = 0;
                for p2 in neighbors_of((x, y)) {
                    if p_live[&p2] {
                        ct += 1;
                    }
                }
                let live = match p_live[&(x, y)] {
                    true => (3 <= ct && ct <= 4),
                    false => (ct == 3),
                };
                if live {
                    let idx = y * mx + x;
                    res |= (1 << idx);
                }
            }

            map[masked] = res;
        }

        debug_log(format!("Adding mask of {} bits covering {} bits", mask_ps.len(), covered.len()));

        (mask, map)
    };

    let mut acc = Vec::new();
    let mut left: Vec<_> = (0..mx).flat_map(|x| (0..my).map(move |y| (x, y))).collect();
    while !left.is_empty() {
        let mut mask_ps = HashSet::new();

        for &p in left.iter() {
            let mut mask_ps2 = mask_ps.clone();
            for p2 in neighbors_of(p) {
                mask_ps2.insert(p2);
            }
            if mask_ps2.len() > max_mask_bits {
                break;
            }
            mask_ps = mask_ps2;
        }

        let mut covered = Vec::new();
        let mut uncovered = Vec::new();

        for p in left {
            if neighbors_of(p).all(|p2| mask_ps.contains(&p2)) {
                covered.push(p);
            }
            else {
                uncovered.push(p);
            }
        }

        assert!(!covered.is_empty());

        acc.push(compute_masks1(mask_ps, covered));

        left = uncovered;
    }
    acc
}

fn main_rand(min_area: isize, max_area: isize) {
    let threads = 8;

    let mut lattices = HashMap::new();
    for mx in 5..=max_area {
        for my in 5..=(max_area / mx) {
            let sz = mx * my;
            if sz < min_area {
                continue;
            }
            let lats = lattices.entry(sz).or_insert_with(|| Vec::new());
            for syx in 0..=(mx / 2) {
                lats.push((mx, my, syx));
            }
        }
    }
    let lattices: Vec<_> = lattices.into_iter().map(|(_, lats)| lats).collect();

    let (tx, rx) = std::sync::mpsc::sync_channel(1024);

    let now = std::time::Instant::now();
    let heartbeats: Vec<_> = (0..threads).map(|_| Arc::new(Mutex::new(now))).collect();

    crossbeam::scope(|sc| {
        for n in 0..threads {
            let lattices = &lattices;
            let tx = tx.clone();
            let heartbeats = &heartbeats;
            sc.spawn(move |_| {
                let mut already = HashSet::new();
                loop {
                    for result in main_rand1(lattices) {
                        if already.contains(&result) {
                            continue;
                        }
                        tx.send(Some(result.clone())).unwrap();
                        already.insert(result);
                    }

                    {
                        let mut hb = heartbeats[n].lock().unwrap();
                        *hb = std::time::Instant::now();
                    }
                }
            });
        }

        {
            let tx = tx.clone();
            sc.spawn(move |_| {
                loop {
                    std::thread::sleep(Duration::from_millis(60000));

                    tx.send(None).unwrap();
                }
            });
        }

        let mut already = HashSet::new();
        loop {
            match rx.recv().unwrap() {
                Some(result) => {
                    if already.contains(&result) {
                        continue;
                    }
                    println!("{}", serde_json::to_string(&result).unwrap());
                    already.insert(result);
                }
                None => {
                    let now = std::time::Instant::now();
                    let mut hbs: Vec<_> = heartbeats.iter().enumerate().map(|(n, hb)| (n, *hb.lock().unwrap())).collect();
                    hbs.sort_by_key(|&(n, hb)| (hb, n));
                    let pieces: Vec<_> = hbs.into_iter().map(|(n, hb)| format!("#{}: {:?}", n, now - hb)).collect();
                    debug_log(format!("Heartbeats: {}", pieces.join(", ")));
                }
            }
        }
    }).unwrap();
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

fn main_rand1(lattices: &[Vec<Vec3>]) -> Vec<(Geometry3, Vec<Vec2>)> {
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
        if let Some(t0) = already.get(&gen0) {
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
