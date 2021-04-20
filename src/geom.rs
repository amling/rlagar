pub type Vec3 = (isize, isize, isize);
pub type Vec2 = (isize, isize);
pub type Vec1 = (isize,);
pub type Geometry2 = (Option<Vec2>, (Option<Vec1>, ()));
pub type Geometry3 = (Option<Vec3>, Geometry2);
