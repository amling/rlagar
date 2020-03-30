pub trait ZModule {
    fn zero() -> Self;
    fn mul(&mut self, q: isize);
    fn submul(&mut self, q: isize, b: &Self);
}

impl ZModule for () {
    fn zero() {
    }

    fn mul(&mut self, _q: isize) {
    }

    fn submul(&mut self, _q: isize, _b: &Self) {
    }
}

impl ZModule for isize {
    fn zero() -> Self {
        0
    }

    fn mul(&mut self, q: isize) {
        *self *= q;
    }

    fn submul(&mut self, q: isize, b: &Self) {
        *self -= q * *b;
    }
}

impl<A: ZModule, B: ZModule> ZModule for (A, B) {
    fn zero() -> Self {
        (A::zero(), B::zero())
    }

    fn mul(&mut self, q: isize) {
        self.0.mul(q);
        self.1.mul(q);
    }

    fn submul(&mut self, q: isize, b: &Self) {
        self.0.submul(q, &b.0);
        self.1.submul(q, &b.1);
    }
}

pub fn egcd<R: ZModule>(mut a: isize, mut b: isize, mut ra: R, mut rb: R) -> (isize, R, R) {
    egcd_mut(&mut a, &mut b, &mut ra, &mut rb);
    (b, ra, rb)
}

fn egcd_mut<R: ZModule>(a: &mut isize, b: &mut isize, ra: &mut R, rb: &mut R) {
    if *a < 0 {
        *a *= -1;
        ra.mul(-1);
    }

    if *b < 0 {
        *b *= -1;
        rb.mul(-1);
    }

    while *a > 0 {
        let q = *b / *a;

        *b -= q * *a;
        rb.submul(q, &ra);

        std::mem::swap(a, b);
        std::mem::swap(ra, rb);
    }
}

pub trait LatticeCanonicalizable: ZModule + Sized {
    type Output;

    fn canonicalize(vs: Vec<Self>) -> Self::Output;
}

impl LatticeCanonicalizable for () {
    type Output = ();

    fn canonicalize(_: Vec<Self>) {
    }
}

impl<S: LatticeCanonicalizable> LatticeCanonicalizable for (S, isize) {
    type Output = (Option<(S, isize)>, S::Output);

    fn canonicalize(mut vs: Vec<(S, isize)>) -> (Option<(S, isize)>, S::Output) {
        let mut l = (S::zero(), 0);

        for v in vs.iter_mut() {
            let (sv, nv) = v;
            let (sl, nl) = &mut l;
            egcd_mut(nv, nl, sv, sl);
        }

        let l = match l.1 != 0 {
            true => Some(l),
            false => None,
        };

        let vs = vs.into_iter().map(|(s, n)| {
            assert_eq!(0, n);
            s
        }).collect();

        // TODO: canonicalize l.0 w.r.t. output of S::canonicalize(vs)

        (l, S::canonicalize(vs))
    }
}
