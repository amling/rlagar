#![allow(unused_parens)]

#[macro_use]
extern crate ars_ds;

use crossbeam::queue::PopError;
use crossbeam::queue::SegQueue;
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

mod ana;
mod flags;
mod geom;
mod lattice;
mod misc;
mod rando;

use flags::Flags;
use flags::SimpleFlags;
use geom::Geometry3;
use geom::Vec2;
use geom::Vec3;
use lattice::CanonicalLattice;
use lattice::LatticeCanonicalizable;

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

    if cmd == "rand" {
        let min_area = args.next().unwrap().parse().unwrap();
        let max_area = args.next().unwrap().parse().unwrap();
        return rando::main_rand(min_area, max_area);
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

fn all_maybe_map<T, R>(i: impl Iterator<Item=T>, mut f: impl FnMut(T) -> Option<R>) -> Option<Vec<R>> {
    let mut ret = Vec::new();
    for t in i {
        if let Some(r) = f(t) {
            ret.push(r);
            continue;
        }
        return None;
    }
    Some(ret)
}

fn gens(mx: isize, my: isize, syx: isize) {
    let masks = compute_masks((mx, my, syx));

    let maybe_masks = all_maybe_map(masks.iter(), |m| {
        match &m[..] {
            [] => panic!(),
            &[m] => Some(m),
            _ => None,
        }
    });
    if let Some(masks) = maybe_masks {
        return gens2(mx, my, syx, &masks);
    }

    let maybe_masks = all_maybe_map(masks.iter(), |m| {
        match &m[..] {
            [] => panic!(),
            &[m] => Some((m, 0)),
            &[m1, m2] => Some((m1, m2)),
            _ => None,
        }
    });
    if let Some(masks) = maybe_masks {
        return gens2(mx, my, syx, &masks);
    }

    let maybe_masks = all_maybe_map(masks.iter(), |m| {
        match &m[..] {
            [] => panic!(),
            &[m] => Some((m, 0, 0)),
            &[m1, m2] => Some((m1, m2, 0)),
            &[m1, m2, m3] => Some((m1, m2, m3)),
            _ => None,
        }
    });
    if let Some(masks) = maybe_masks {
        return gens2(mx, my, syx, &masks);
    }

    gens2(mx, my, syx, &masks);
}

trait Mask: Send + Sync {
    fn count(&self, s: u64) -> u32;
}

impl Mask for u64 {
    fn count(&self, s: u64) -> u32 {
        (s & self).count_ones()
    }
}

impl Mask for (u64, u64) {
    fn count(&self, s: u64) -> u32 {
        (s & self.0).count_ones() + (s & self.1).count_ones()
    }
}

impl Mask for (u64, u64, u64) {
    fn count(&self, s: u64) -> u32 {
        (s & self.0).count_ones() + (s & self.1).count_ones() + (s & self.2).count_ones()
    }
}

impl Mask for Vec<u64> {
    fn count(&self, s: u64) -> u32 {
        let mut ct = 0;
        for mask in self {
            ct += (s & mask).count_ones();
        }
        ct
    }
}

fn gens2(mx: isize, my: isize, syx: isize, masks: &Vec<impl Mask>) {
    let lattice = (mx, my, syx);
    misc::debug_time(format!("lattice {:?}", lattice), || {
        let n = mx * my;

        let threads = 8;
        let workunit_bits = 20.min(n);

        let geometry2 = (Some((syx, my)), (Some((mx,)), ()));

        let flags = SimpleFlags::new(1 << n);
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
                                search(masks, flags, s0, results);
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
                            misc::debug_log(format!("Completed {}/{}, estimated {:?} left...", (total - now.1), total, left));
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

                s1 = tick(masks, s1);
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
            ana::ana2(mx, my, syx, period, gen0)
        });

        let results: BTreeSet<_> = results.collect();

        misc::debug_log(format!("Lattice {:?} => {} results", lattice, results.len()));

        for res in results {
            println!("{}", serde_json::to_string(&res).unwrap());
        }
    });
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

fn search(masks: &Vec<impl Mask>, flags: &impl Flags, s0: u64, results: &mut HashSet<(u64, isize)>) {
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

fn tick(masks: &Vec<impl Mask>, s0: u64) -> u64 {
    let mut s1 = 0;
    for (idx, mask) in masks.iter().enumerate() {
        let ct = mask.count(s0);
        let self_ct = ((s0 >> idx) as u32) & 1;
        let magic_ct = ct * 2 + self_ct;
        if 5 <= magic_ct && magic_ct <= 7 {
            s1 |= (1 << idx);
        }
    }
    s1
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
