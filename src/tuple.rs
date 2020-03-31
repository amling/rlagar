pub trait TupleEnd<B> {
    type F;

    fn split_tuple_end(zelf: Self) -> (Self::F, B);
    fn join_tuple_end(f: Self::F, b: B) -> Self;
}

impl<A> TupleEnd<A> for A {
    type F = ();

    fn split_tuple_end(a: Self) -> ((), A) {
        ((), a)
    }

    fn join_tuple_end(_: (), a: A) -> A {
        a
    }
}

impl<A, B> TupleEnd<B> for (A, B) {
    type F = A;

    fn split_tuple_end((a, b): (A, B)) -> (A, B) {
        (a, b)
    }

    fn join_tuple_end(a: A, b: B) -> (A, B) {
        (a, b)
    }
}

impl<A, B, C> TupleEnd<C> for (A, B, C) {
    type F = (A, B);

    fn split_tuple_end((a, b, c): (A, B, C)) -> ((A, B), C) {
        ((a, b), c)
    }

    fn join_tuple_end((a, b): (A, B), c: C) -> (A, B, C) {
        (a, b, c)
    }
}
