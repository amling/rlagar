mod bits;

fn egcd(a: isize, b: isize) -> (isize, isize, isize) {
    let mut a = a.abs();
    let mut b = b.abs();

    let mut sa = 1;
    let mut ta = 0;

    let mut sb = 0;
    let mut tb = 1;

    while a > 0 {
        let q = b / a;

        sb -= sa * q;
        tb -= ta * q;
        b -= a * q;

        std::mem::swap(&mut sa, &mut sb);
        std::mem::swap(&mut ta, &mut tb);
        std::mem::swap(&mut a, &mut b);
    }

    (b, sb, tb)
}

fn main() {
    let n = 12;

    let mut lattices = Vec::new();
    for w in 1..=n {
        if n % w != 0 {
            continue;
        }

        let h = n / w;

        for sx in 0..w {
            if 2 * sx > w {
                // don't bother with shifts over half, instead we'll find whatever flipped
                continue;
            }

            let (d, s, t) = egcd(sx, w);
            // d = s * sx + t * w

            // (w, 0), (sx, h)
            // (w, 0), (sx, h), (sx * w / d, h * w / d), (s * sx + t * w, s * h + t * 0)
            // (w, 0), (sx, h), (sx/d * w, h * w / d), (d, s * h)
            // (w, 0), (sx, h), (0, h * w / d), (d, s * h)

            // the order stepping vertically, what w would be if we transposed
            let wt = h * w / d;
            // the rest...
            let ht = d;
            // what sx would be if we transposed
            let sxt = s * t;
            // grumble grumble stupid mod
            let sxt = ((sxt % wt) + wt) % wt;
            // if smaller flipped...
            let sxt = sxt.min(wt - sxt);

            if (-wt, ht, sxt) < (-w, ht, sx) {
                continue;
            }

            lattices.push((w, h, sx));
        }
    }

    println!("Lattices {:?}", lattices);
}
