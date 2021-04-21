use core::hash::Hash;
use rand::seq::SliceRandom;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Duration;

use crate::ana;
use crate::engine;
use crate::geom;
use crate::misc;

use engine::Engine;
use engine::MaskEngine;
use engine::SetEngine;
use geom::Geometry3;
use geom::Vec2;
use geom::Vec3;

enum RandEngine {
    MaskU64(MaskEngine<u64>),
    MaskU128(MaskEngine<u128>),
    Set(SetEngine),
}

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
    let lattices = lattices.into_iter().map(|(_, lats)| lats);
    let lattices: Vec<_> = lattices.map(|lats| {
        lats.into_iter().map(|(mx, my, syx)| {
            let engine = if let Some(engine) = MaskEngine::<Vec<u64>>::compile(mx, my, syx).and_then(|e| e.remask_vec_single()) {
                RandEngine::MaskU64(engine)
            }
            else if let Some(engine) = MaskEngine::<Vec<u128>>::compile(mx, my, syx).and_then(|e| e.remask_vec_single()) {
                RandEngine::MaskU128(engine)
            }
            else {
                RandEngine::Set(SetEngine::compile(mx, my, syx))
            };

            ((mx, my, syx), engine)
        }).collect()
    }).collect();

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
                    main_rand1(lattices, |result| {
                        if !already.contains(&result) {
                            tx.send(Some(result.clone())).unwrap();
                            already.insert(result);
                        }

                        let mut hb = heartbeats[n].lock().unwrap();
                        *hb = std::time::Instant::now();
                    });
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

fn main_rand1(lattices: &[Vec<(Vec3, RandEngine)>], f: impl FnMut((Geometry3, Vec<Vec2>))) {
    let lats = lattices.choose(&mut rand::thread_rng()).unwrap();
    let &((mx, my, syx), ref engine) = lats.choose(&mut rand::thread_rng()).unwrap();

    match engine {
        RandEngine::MaskU64(engine) => main_rand2(mx, my, syx, engine, f),
        RandEngine::MaskU128(engine) => main_rand2(mx, my, syx, engine, f),
        RandEngine::Set(engine) => main_rand2(mx, my, syx, engine, f),
    }
}

trait RandHashable {
    type H: Hash + Eq;

    fn to_hashable(&self) -> Self::H;
}

impl RandHashable for u64 {
    type H = u64;

    fn to_hashable(&self) -> u64 {
        *self
    }
}

impl RandHashable for u128 {
    type H = u128;

    fn to_hashable(&self) -> u128 {
        *self
    }
}

impl RandHashable for HashSet<Vec2> {
    type H = Vec<Vec2>;

    fn to_hashable(&self) -> Vec<Vec2> {
        let mut ret: Vec<_> = self.iter().cloned().collect();
        ret.sort();
        ret
    }
}

fn main_rand2<S: RandHashable>(mx: isize, my: isize, syx: isize, engine: &impl Engine<S>, mut f: impl FnMut((Geometry3, Vec<Vec2>))) {
    let mut s0 = engine.rand();
    let mut t = 0isize;
    let mut already = HashMap::new();
    loop {
        let h0 = s0.to_hashable();
        if let Some(t0) = already.get(&h0) {
            for r in ana::ana2(mx, my, syx, t - t0, engine.decode(&s0)) {
                f(r)
            }
            return;
        }
        already.insert(h0, t);

        s0 = engine.tick(&s0);
        t += 1;
    }
}
