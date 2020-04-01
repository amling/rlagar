use ars_aa::zmodule::ZModule;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;

pub fn search_lattice_links<N: Hash + Eq + Clone, R: ZModule + Clone>(links: &HashMap<N, HashSet<(N, R)>>, n: N) -> (HashSet<N>, Vec<R>) {
    let mut found = HashMap::new();
    find_connected(links, &mut found, n, R::zero());

    let mut acc = Vec::new();
    for (n1, links2) in links {
        for (n2, r) in links2 {
            if let Some(r1) = found.get(n1) {
                if let Some(r2) = found.get(n2) {
                    let mut label = r1.clone();
                    label.addmul(1, &r);
                    label.addmul(-1, &r2);
                    acc.push(label);
                }
            }
        }
    }

    let found = found.into_iter().map(|(n, _)| n).collect();

    (found, acc)
}

fn find_connected<N: Hash + Eq + Clone, R: ZModule + Clone>(links: &HashMap<N, HashSet<(N, R)>>, found: &mut HashMap<N, R>, n: N, r: R) {
    if let Some(_) = found.get(&n) {
        return;
    }

    found.insert(n.clone(), r.clone());

    for (n2, rd) in links.get(&n).unwrap() {
        let mut r2 = r.clone();
        r2.addmul(1, rd);
        find_connected(links, found, n2.clone(), r2);
    }
}
