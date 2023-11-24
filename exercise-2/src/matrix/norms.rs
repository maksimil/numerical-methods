use crate::{
    basic::{Numerical, OtherNumericalOps},
    representation::repr_ref,
};

use super::{
    column::{ColumnOf, ColumnRef},
    traits::MatrixRef,
    transpose::MatrixTranspose,
};

pub trait NormedColumn: ColumnRef {
    fn norm_one(&self) -> Self::Scalar;
    fn norm_inf(&self) -> Self::Scalar;
}

impl<T> NormedColumn for T
where
    T: ColumnRef,
    T::Scalar: Numerical,
{
    fn norm_one(&self) -> Self::Scalar
    where
        Self::Scalar: Numerical,
    {
        (0..self.dimension()).map(|i| self.at(i).abs_trait()).sum()
    }

    fn norm_inf(&self) -> Self::Scalar
    where
        Self::Scalar: Numerical,
    {
        (0..self.dimension())
            .map(|i| self.at(i).abs_trait())
            .fold(Self::Scalar::zero(), |v, m| if v > m { v } else { m })
    }
}

pub trait NormedMatrix: MatrixRef {
    fn norm_one(&self) -> Self::Scalar;
    fn norm_inf(&self) -> Self::Scalar;
}

impl<T> NormedMatrix for T
where
    T: MatrixRef,
    T::Scalar: Numerical,
{
    fn norm_one(&self) -> Self::Scalar {
        (0..self.dimension())
            .map(|j| ColumnOf::new(repr_ref::<T>(&self), j).norm_one())
            .fold(Self::Scalar::zero(), |v, m| if v > m { v } else { m })
    }

    fn norm_inf(&self) -> Self::Scalar {
        MatrixTranspose::from(repr_ref::<T>(&self)).norm_one()
    }
}
