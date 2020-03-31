use crate::tuple;

use tuple::TupleEnd;

pub trait ZModule: Eq + Clone {
    fn zero() -> Self;
    fn mul(&mut self, q: isize);
    fn addmul(&mut self, q: isize, b: &Self);
}

impl ZModule for () {
    fn zero() {
    }

    fn mul(&mut self, _q: isize) {
    }

    fn addmul(&mut self, _q: isize, _b: &Self) {
    }
}

impl ZModule for isize {
    fn zero() -> Self {
        0
    }

    fn mul(&mut self, q: isize) {
        *self *= q;
    }

    fn addmul(&mut self, q: isize, b: &Self) {
        *self += q * *b;
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

    fn addmul(&mut self, q: isize, b: &Self) {
        self.0.addmul(q, &b.0);
        self.1.addmul(q, &b.1);
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
        rb.addmul(-q, &ra);

        std::mem::swap(a, b);
        std::mem::swap(ra, rb);
    }
}

pub trait Canonicalizes<S: ZModule> {
    fn canonicalize(&self, s: S) -> S;

    fn canonicalize_delta(&self, s: S) -> (S, S) {
        let mut s2 = s.clone();
        s2 = self.canonicalize(s2);
        let mut sd = s2.clone();
        sd.addmul(-1, &s);
        (s, sd)
    }
}

pub trait LatticeCanonicalizable: ZModule + Sized {
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

impl<S: LatticeCanonicalizable, T: ZModule + TupleEnd<isize, F=S>> LatticeCanonicalizable for T {
    type Output = (Option<T>, S::Output);

    fn canonicalize(vs: Vec<T>) -> (Option<T>, S::Output) {
        let mut l = (S::zero(), 0);

        let mut vs: Vec<_> = vs.into_iter().map(|t| T::split_tuple_end(t)).collect();

        for v in vs.iter_mut() {
            let (sv, nv) = v;
            let (sl, nl) = &mut l;
            egcd_mut(nv, nl, sv, sl);
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

impl<S: LatticeCanonicalizable, T: ZModule + TupleEnd<isize, F=S>> Canonicalizes<T> for (Option<T>, S::Output) {
    fn canonicalize(&self, t: T) -> T {
        let (mut s, mut n) = T::split_tuple_end(t);
        if let Some(t) = &self.0 {
            let (s1, n1) = T::split_tuple_end(t.clone());

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

type Vec0 = ();
type Vec1 = (Vec0, isize);
type Vec2 = (Vec1, isize);
type Vec3 = (Vec2, isize);

type Geom0 = ();
type Geom1 = (Option<Vec1>, Geom0);
type Geom2 = (Option<Vec2>, Geom1);
type Geom3 = (Option<Vec3>, Geom2);

pub fn vec0() -> Vec0 {
    ()
}

pub fn vec1(x: isize) -> Vec1 {
    (vec0(), x)
}

pub fn vec2(x: isize, y: isize) -> Vec2 {
    (vec1(x), y)
}

pub fn vec3(x: isize, y: isize, t: isize) -> Vec3 {
    (vec2(x, y), t)
}

pub fn unvec1(v: Vec1) -> isize {
    let ((), x) = v;
    x
}

pub fn unvec2(v: Vec2) -> (isize, isize) {
    let (((), x), y) = v;
    (x, y)
}

pub fn unvec3(v: Vec3) -> (isize, isize, isize) {
    let ((((), x), y), t) = v;
    (x, y, t)
}

pub fn geom1(mx: isize) -> Geom1 {
    (Some(vec1(mx)), ())
}

pub fn geom2(mx: isize, my: isize, syx: isize) -> Geom2 {
    (Some(vec2(syx, my)), geom1(mx))
}

pub fn geom3(mx: isize, my: isize, syx: isize, mt: isize, stx: isize, sty: isize) -> Geom3 {
    (Some(vec3(stx, sty, mt)), geom2(mx, my, syx))
}
