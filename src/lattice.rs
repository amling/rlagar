pub trait ZModule {
    fn mul(&mut self, q: isize);
    fn submul(&mut self, q: isize, b: &Self);
}

impl ZModule for () {
    fn mul(&mut self, _q: isize) {
    }

    fn submul(&mut self, _q: isize, _b: &Self) {
    }
}

impl ZModule for isize {
    fn mul(&mut self, q: isize) {
        *self *= q;
    }

    fn submul(&mut self, q: isize, b: &Self) {
        *self -= q * *b;
    }
}

impl<A: ZModule, B: ZModule> ZModule for (A, B) {
    fn mul(&mut self, q: isize) {
        self.0.mul(q);
        self.1.mul(q);
    }

    fn submul(&mut self, q: isize, b: &Self) {
        self.0.submul(q, &b.0);
        self.1.submul(q, &b.1);
    }
}

pub fn egcd<R: ZModule>(a: isize, b: isize, ra: R, rb: R) -> (isize, R, R) {
    let mut a = a;
    let mut ra = ra;
    if a < 0 {
        a = -a;
        ra.mul(-1);
    }

    let mut b = b;
    let mut rb = rb;
    if b < 0 {
        b = -b;
        rb.mul(-1);
    }

    while a > 0 {
        let q = b / a;

        b -= q * a;
        rb.submul(q, &ra);

        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut ra, &mut rb);
    }

    (b, ra, rb)
}
