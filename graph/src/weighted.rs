use ars_aa::zmodule::ZModule;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::hash::Hash;

// Finds generators for the space of cycle weights reachable from the given node.
pub fn find_cycle_generators<N: Hash + Clone + Eq, R: ZModule + Clone>(links: &HashMap<N, HashSet<(N, R)>>, n: N) -> Vec<R> {
    let connected = find_connected(links, n);

    let mut acc = Vec::new();
    for (n1, links2) in links {
        for (n2, r) in links2 {
            if let Some(r1) = connected.get(n1) {
                if let Some(r2) = connected.get(n2) {
                    let mut label = r1.clone();
                    label.addmul(1, r);
                    label.addmul(-1, r2);
                    acc.push(label);
                }
            }
        }
    }
    acc
}

// Finds the set of connected nodes and the weight of some arbitrary path to them
pub fn find_connected<N: Hash + Clone + Eq, R: ZModule + Clone>(links: &HashMap<N, HashSet<(N, R)>>, n: N) -> HashMap<N, R> {
    let mut acc = HashMap::new();
    let mut q = VecDeque::new();
    q.push_back((n, R::zero()));
    loop {
        let (n, r) = match q.pop_front() {
            Some(e) => e,
            None => {
                return acc;
            }
        };
        if let Some(_) = acc.get(&n) {
            continue;
        }

        acc.insert(n.clone(), r.clone());

        for (n2, rd) in links.get(&n).unwrap() {
            let mut r2 = r.clone();
            r2.addmul(1, rd);
            q.push_back((n2.clone(), r2));
        }
    }
}
