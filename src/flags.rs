use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering;

pub struct SimpleFlags {
    vec: Vec<AtomicU64>,
}

impl SimpleFlags {
    pub fn new(sz: u64) -> Self {
        let sz = (sz + 63) / 64;
        SimpleFlags {
            vec: (0..sz).map(|_| AtomicU64::new(0)).collect(),
        }
    }

    pub fn get(&self, idx: u64) -> bool {
        let shift = idx % 64;
        let idx = (idx / 64) as usize;
        let mask = 1 << shift;
        self.vec[idx].load(Ordering::Relaxed) & mask == mask
    }

    pub fn set(&self, idx: u64) {
        let shift = idx % 64;
        let idx = (idx / 64) as usize;
        let mask = 1 << shift;

        loop {
            let v = self.vec[idx].load(Ordering::Relaxed);
            if self.vec[idx].compare_and_swap(v, v | mask, Ordering::Relaxed) == v {
                return;
            }
        }
    }
}
