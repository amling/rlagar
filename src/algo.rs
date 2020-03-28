use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;

pub trait LatticeResult: Clone {
    fn zero() -> Self;
    fn neg(self) -> Self;
    fn add(self, rhs: Self) -> Self;
}

impl LatticeResult for isize {
    fn zero() -> Self {
        0
    }

    fn neg(self) -> Self {
        -self
    }

    fn add(self, rhs: Self) -> Self {
        self + rhs
    }
}

impl<A: LatticeResult, B: LatticeResult> LatticeResult for (A, B) {
    fn zero() -> Self {
        (A::zero(), B::zero())
    }

    fn neg(self) -> Self {
        (self.0.neg(), self.1.neg())
    }

    fn add(self, rhs: Self) -> Self {
        (self.0.add(rhs.0), self.1.add(rhs.1))
    }
}

impl<A: LatticeResult, B: LatticeResult, C: LatticeResult> LatticeResult for (A, B, C) {
    fn zero() -> Self {
        (A::zero(), B::zero(), C::zero())
    }

    fn neg(self) -> Self {
        (self.0.neg(), self.1.neg(), self.2.neg())
    }

    fn add(self, rhs: Self) -> Self {
        (self.0.add(rhs.0), self.1.add(rhs.1), self.2.add(rhs.2))
    }
}

pub fn search_lattice_links<N: Hash + Eq + Clone, R: LatticeResult>(links: &HashMap<N, HashSet<(N, R)>>, n: N) -> (HashSet<N>, Vec<R>) {
    let mut found = HashMap::new();
    find_connected(links, &mut found, n, R::zero());

    let mut acc = Vec::new();
    for (n1, links2) in links {
        for (n2, r) in links2 {
            if let Some(r1) = found.get(n1) {
                if let Some(r2) = found.get(n2) {
                    acc.push(r1.clone().add(r.clone()).add(r2.clone().neg()));
                }
            }
        }
    }

    let found = found.into_iter().map(|(n, _)| n).collect();

    (found, acc)
}

fn find_connected<N: Hash + Eq + Clone, R: LatticeResult>(links: &HashMap<N, HashSet<(N, R)>>, found: &mut HashMap<N, R>, n: N, r: R) {
    if let Some(_) = found.get(&n) {
        return;
    }

    found.insert(n.clone(), r.clone());

    for (n2, r2) in links.get(&n).unwrap() {
        find_connected(links, found, n2.clone(), r.clone().add(r2.clone()));
    }
}
