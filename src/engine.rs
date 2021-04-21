use ars_aa::lattice::LatticeCanonicalizer;
use rand::Rng;
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
    fn rand(bits: usize) -> Self;
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

            fn rand(bits: usize) -> Self {
                assert!(bits <= Self::sz());
                let mut ret = rand::thread_rng().gen();
                if bits < Self::sz() {
                    ret &= ((1 as Self) << (bits) - 1);
                }
                ret
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

#[derive(Clone)]
#[derive(Copy)]
struct LatticeShape {
    mx: isize,
    my: isize,
    syx: isize,
}

impl LatticeShape {
    pub fn new(mx: isize, my: isize, syx: isize) -> Self {
        LatticeShape {
            mx: mx,
            my: my,
            syx: syx,
        }
    }

    pub fn geom2(&self) -> Geometry2 {
        (Some((self.syx, self.my)), (Some((self.mx,)), ()))
    }
}

pub trait Engine<S> {
    fn tick(&self, s0: &S) -> S;
    fn rand(&self) -> S;
    fn decode(&self, s: &S) -> HashSet<Vec2>;
}

pub struct MaskEngine<M>(LatticeShape, Vec<M>);

impl<S: BitState> MaskEngine<Vec<S>> {
    pub fn compile(mx: isize, my: isize, syx: isize) -> Option<Self> {
        if mx * my > (S::sz() as isize) {
            return None;
        }
        let shape = LatticeShape::new(mx, my, syx);
        let geometry2 = shape.geom2();

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
        Some(MaskEngine(shape, acc))
    }
}

impl<M> MaskEngine<M> {
    pub fn remask<M2>(&self, mut f: impl FnMut(&M) -> Option<M2>) -> Option<MaskEngine<M2>> {
        let mut ret = Vec::new();
        for m in self.1.iter() {
            match f(m) {
                Some(m2) => {
                    ret.push(m2);
                }
                None => {
                    return None
                }
            }
        }
        Some(MaskEngine(self.0, ret))
    }
}

impl<S: BitState + Copy> MaskEngine<Vec<S>> {
    pub fn remask_vec_single(&self) -> Option<MaskEngine<S>> {
        self.remask(|m| {
            match &m[..] {
                [] => panic!(),
                &[m] => Some(m),
                _ => None,
            }
        })
    }
}

impl<S: BitState, M: Mask<S>> Engine<S> for MaskEngine<M> {
    fn tick(&self, s0: &S) -> S {
        let mut s1 = S::zero();
        for (idx, mask) in self.1.iter().enumerate() {
            let ct = mask.count(&s0);
            let self_ct = s0.read_bit_u32(idx);
            let magic_ct = ct * 2 + self_ct;
            if 5 <= magic_ct && magic_ct <= 7 {
                s1.set_bit_true(idx);
            }
        }
        s1
    }

    fn rand(&self) -> S {
        S::rand((self.0.mx * self.0.my) as usize)
    }

    fn decode(&self, s: &S) -> HashSet<Vec2> {
        let mut ret = HashSet::new();
        for x in 0..self.0.mx {
            for y in 0..self.0.my {
                let idx = y * self.0.mx + x;
                if s.read_bit_bool(idx as usize) {
                    ret.insert((x, y));
                }
            }
        }
        ret
    }
}

pub struct SetEngine(LatticeShape, Geometry2);

impl SetEngine {
    pub fn compile(mx: isize, my: isize, syx: isize) -> Self {
        let shape = LatticeShape::new(mx, my, syx);
        SetEngine(shape, shape.geom2())
    }
}

impl Engine<HashSet<Vec2>> for SetEngine {
    fn tick(&self, s0: &HashSet<Vec2>) -> HashSet<Vec2> {
        let mut cts = HashMap::new();
        for &(x, y) in s0 {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    let x2 = x + dx;
                    let y2 = y + dy;
                    let (cx2, cy2) = self.1.canonicalize((x2, y2));
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

    fn rand(&self) -> HashSet<Vec2> {
        let mut s = HashSet::new();
        for x in 0..self.0.mx {
            for y in 0..self.0.my {
                if rand::thread_rng().gen() {
                    s.insert((x, y));
                }
            }
        }
        s
    }

    fn decode(&self, s: &HashSet<Vec2>) -> HashSet<Vec2> {
        s.clone()
    }
}
