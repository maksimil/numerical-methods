use std::{
    iter::Sum,
    ops::{Add, Div, Mul, Neg, Sub},
};

pub type Index = usize;
pub const INDEX_NOT_FOUND: Index = usize::MAX;

pub trait OtherNumericalOps {
    fn abs_trait(&self) -> Self;
    fn zero() -> Self;
}

macro_rules! impl_other_numerical_ops {
    ($ty:ident) => {
        impl OtherNumericalOps for $ty {
            fn abs_trait(&self) -> Self {
                self.abs()
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
        + Mul<Self, Output = Self>
        + Div<Self, Output = Self>
        + OtherNumericalOps,
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
        + Mul<Self, Output = Self>
        + Div<Self, Output = Self>
        + OtherNumericalOps
{
}
