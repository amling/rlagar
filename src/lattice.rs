use ars_aa::zmodule::ZModule;
use ars_ds::tuple::CTupleEnd;

is_tuple_trait!(IsTuple);

pub trait Canonicalizes<S: ZModule + Clone> {
    fn canonicalize(&self, s: S) -> S;

    fn canonicalize_delta(&self, s: S) -> (S, S) {
        let mut s2 = s.clone();
        s2 = self.canonicalize(s2);
        let mut sd = s2.clone();
        sd.addmul(-1, &s);
        (s, sd)
    }
}

pub trait LatticeCanonicalizable: ZModule + Clone + Sized {
    type Output: Canonicalizes<Self>;

    fn canonicalize(vs: Vec<Self>) -> Self::Output;
}

impl LatticeCanonicalizable for () {
    type Output = ();

    fn canonicalize(_: Vec<Self>) {
    }
}

impl Canonicalizes<()> for () {
    fn canonicalize(&self, _: ()) -> () {
    }
}

impl<S: LatticeCanonicalizable + Eq, T: ZModule + CTupleEnd<F=S, B=isize> + IsTuple + Clone> LatticeCanonicalizable for T {
    type Output = (Option<T>, S::Output);

    fn canonicalize(vs: Vec<T>) -> (Option<T>, S::Output) {
        let mut l = (S::zero(), 0);

        let mut vs: Vec<_> = vs.into_iter().map(|t| T::split_tuple_end(t)).collect();

        for v in vs.iter_mut() {
            let (sv, nv) = v;
            let (sl, nl) = &mut l;
            ars_aa::misc::egcd_mut(nv, nl, sv, sl);
        }

        let vs = vs.into_iter().map(|(s, n)| {
            assert_eq!(0, n);
            s
        }).collect();

        let so = S::canonicalize(vs);

        let l = match l.1 != 0 {
            true => Some(T::join_tuple_end(so.canonicalize(l.0), l.1)),
            false => {
                assert!(S::zero() == l.0);
                None
            },
        };

        (l, so)
    }
}

impl<S: LatticeCanonicalizable, T: ZModule + CTupleEnd<F=S, B=isize> + Clone> Canonicalizes<T> for (Option<T>, S::Output) {
    fn canonicalize(&self, t: T) -> T {
        let (mut s, mut n) = T::split_tuple_end(t);
        if let Some(t) = &self.0 {
            let (s1, n1) = T::split_tuple_end(t.clone());

            // unbelievably, this absolutely smokes using div_euclid
            while n < 0 {
                n += n1;
                s.addmul(1, &s1);
            }
            while n >= n1 {
                n -= n1;
                s.addmul(-1, &s1);
            }
        }
        s = self.1.canonicalize(s);
        T::join_tuple_end(s, n)
    }
}
