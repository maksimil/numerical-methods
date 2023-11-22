use std::{
    iter::Sum,
    ops::{Add, Div, DivAssign, Mul, Neg, Sub, SubAssign},
};

pub type Index = usize;
pub const INDEX_NOT_FOUND: Index = usize::MAX;

pub trait OtherNumericalOps {
    fn abs_trait(&self) -> Self;
    fn sqrt_trait(&self) -> Self;
    fn zero() -> Self;
}

macro_rules! impl_other_numerical_ops {
    ($ty:ident) => {
        impl OtherNumericalOps for $ty {
            fn abs_trait(&self) -> Self {
                self.abs()
            }
            fn sqrt_trait(&self) -> Self {
                self.sqrt()
            }
            fn zero() -> Self {
                0.0
            }
        }
    };
}

impl_other_numerical_ops!(f32);
impl_other_numerical_ops!(f64);

pub trait Numerical
where
    Self: Default
        + Clone
        + Copy
        + PartialOrd<Self>
        + Neg<Output = Self>
        + Add<Self, Output = Self>
        + Sum
        + Sub<Self, Output = Self>
        + SubAssign<Self>
        + Mul<Self, Output = Self>
        + Div<Self, Output = Self>
        + DivAssign<Self>
        + OtherNumericalOps
        + From<i32>,
{
}

impl<T> Numerical for T where
    Self: Default
        + Clone
        + Copy
        + PartialOrd<Self>
        + Neg<Output = Self>
        + Add<Self, Output = Self>
        + Sum
        + Sub<Self, Output = Self>
        + SubAssign<Self>
        + Mul<Self, Output = Self>
        + Div<Self, Output = Self>
        + DivAssign<Self>
        + OtherNumericalOps
        + From<i32>
{
}
