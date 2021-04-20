use ars_aa::lattice::LatticeCanonicalizable;
use ars_aa::lattice::LatticeCanonicalizer;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::geom;

use geom::Geometry3;
use geom::Vec2;
use geom::Vec3;

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

pub fn ana2(mx: isize, my: isize, syx: isize, ttp: isize, gen0: HashSet<Vec2>) -> Vec<(Geometry3, Vec<Vec2>)> {
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
