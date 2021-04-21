use ars_aa::lattice::LatticeCanonicalizer;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::geom;

use geom::Geometry2;
use geom::Vec2;

pub trait BitState {
    fn zero() -> Self;
    fn sz() -> usize;
    fn read_bit_u32(&self, idx: usize) -> u32;
    fn read_bit_bool(&self, idx: usize) -> bool;
    fn set_bit_true(&mut self, idx: usize);
}

macro_rules! impl_bit_state {
    ($t:ty, $n:expr) => {
        impl BitState for $t {
            fn zero() -> Self {
                0
            }

            fn sz() -> usize {
                $n
            }

            fn read_bit_u32(&self, idx: usize) -> u32 {
                ((self >> idx) & 1) as u32
            }

            fn read_bit_bool(&self, idx: usize) -> bool {
                ((self >> idx) & 1) != 0
            }

            fn set_bit_true(&mut self, idx: usize) {
                *self |= (1 << idx);
            }
        }
    }
}

impl_bit_state!(u64, 64);
impl_bit_state!(u128, 128);

pub trait Mask<S> {
    fn count(&self, s: &S) -> u32;
}

macro_rules! impl_masks {
    ($t:ty) => {
        impl Mask<$t> for $t {
            fn count(&self, s: &$t) -> u32 {
                (s & self).count_ones()
            }
        }

        impl Mask<$t> for ($t, $t) {
            fn count(&self, s: &$t) -> u32 {
                (s & self.0).count_ones() + (s & self.1).count_ones()
            }
        }

        impl Mask<$t> for ($t, $t, $t) {
            fn count(&self, s: &$t) -> u32 {
                (s & self.0).count_ones() + (s & self.1).count_ones() + (s & self.2).count_ones()
            }
        }

        impl Mask<$t> for Vec<$t> {
            fn count(&self, s: &$t) -> u32 {
                let mut ct = 0;
                for mask in self {
                    ct += (s & mask).count_ones();
                }
                ct
            }
        }
    }
}

impl_masks!(u64);
impl_masks!(u128);

pub trait Engine<S> {
    fn tick(&self, s0: &S) -> S;
}

pub struct MaskEngine<M>(Vec<M>);

impl<S: BitState> MaskEngine<Vec<S>> {
    pub fn compile(mx: isize, my: isize, syx: isize) -> Option<Self> {
        if mx * my > (S::sz() as isize) {
            return None;
        }
        let geometry2 = (Some((syx, my)), (Some((mx,)), ()));

        let mut acc = Vec::new();
        for idx in 0..(mx * my) {
            let x = idx % mx;
            let y = idx / mx;
            let mut masks = Vec::<S>::new();
            let mut add_mask = |idx2| {
                for mask in masks.iter_mut() {
                    if !mask.read_bit_bool(idx2) {
                        mask.set_bit_true(idx2);
                        return;
                    }
                }
                masks.push(S::zero());
                masks.last_mut().unwrap().set_bit_true(idx2);
            };
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if (dx, dy) == (0, 0) {
                        continue;
                    }
                    let (x2, y2) = geometry2.canonicalize((x + dx, y + dy));
                    let idx2 = y2 * mx + x2;
                    add_mask(idx2 as usize);
                }
            }
            acc.push(masks);
        }
        Some(MaskEngine(acc))
    }
}

impl<M> MaskEngine<M> {
    pub fn remask<M2>(&self, mut f: impl FnMut(&M) -> Option<M2>) -> Option<MaskEngine<M2>> {
        let mut ret = Vec::new();
        for m in self.0.iter() {
            match f(m) {
                Some(m2) => {
                    ret.push(m2);
                }
                None => {
                    return None
                }
            }
        }
        Some(MaskEngine(ret))
    }
}

impl<S: BitState, M: Mask<S>> Engine<S> for MaskEngine<M> {
    fn tick(&self, s0: &S) -> S {
        let mut s1 = S::zero();
        for (idx, mask) in self.0.iter().enumerate() {
            let ct = mask.count(&s0);
            let self_ct = s0.read_bit_u32(idx);
            let magic_ct = ct * 2 + self_ct;
            if 5 <= magic_ct && magic_ct <= 7 {
                s1.set_bit_true(idx);
            }
        }
        s1
    }
}

struct SetEngine(Geometry2);

impl Engine<HashSet<Vec2>> for SetEngine {
    fn tick(&self, s0: &HashSet<Vec2>) -> HashSet<Vec2> {
        let mut cts = HashMap::new();
        for &(x, y) in s0 {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    let x2 = x + dx;
                    let y2 = y + dy;
                    let (cx2, cy2) = self.0.canonicalize((x2, y2));
                    *cts.entry((cx2, cy2)).or_insert_with(|| 0) += 1;
                }
            }
        }
        let mut s1 = HashSet::new();
        for ((x, y), ct) in cts.into_iter() {
            let living_curr = s0.contains(&(x, y));

            let living_next = match living_curr {
                true => (3 <= ct && ct <= 4),
                false => ct == 3,
            };
            if living_next {
                s1.insert((x, y));
            }
        }

        s1
    }
}
