use rand::Rng;
use rand::seq::SliceRandom;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Duration;

use crate::ana;
use crate::geom;
use crate::lattice;
use crate::misc;

use geom::Geometry3;
use geom::Vec2;
use geom::Vec3;
use lattice::CanonicalLattice;

pub fn main_rand(min_area: isize, max_area: isize) {
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
                    misc::debug_log(format!("Heartbeats: {}", pieces.join(", ")));
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
            return ana::ana2(mx, my, syx, t - t0, gen0.iter().cloned().collect());
        }
        already.insert(gen0.clone(), t);

        let gen1 = compute_step(mx, my, syx, &gen0.iter().cloned().collect());
        let mut gen1: Vec<_> = gen1.into_iter().collect();
        gen1.sort();

        t += 1;
        gen0 = gen1;
    }
}
