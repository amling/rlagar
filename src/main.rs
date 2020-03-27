extern crate crossbeam;

use crossbeam::queue::PopError;
use crossbeam::queue::SegQueue;
use std::sync::atomic::AtomicU64;

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
        let flags: Vec<_> = (0..(1 << (n - 6))).map(|_| AtomicU64::new(0)).collect();
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

        let results: Vec<_> = results.into_iter().flatten().collect();
        eprintln!("{} results", results.len());
    }
}

fn search(lattice: (isize, isize, isize), flags: &Vec<AtomicU64>, s0: u64, results: &mut Vec<(u64, usize)>) {
}
