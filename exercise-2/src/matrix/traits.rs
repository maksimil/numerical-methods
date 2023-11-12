use std::borrow::{Borrow, BorrowMut};

use crate::{basic::Index, representation::Representation};

pub trait MatrixRef<Scalar> {
    fn dimension(&self) -> Index;
    fn at(&self, row: Index, column: Index) -> &Scalar;
}

pub trait MatrixMutRef<Scalar>: MatrixRef<Scalar> {
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar;
}

pub trait MatrixFuncInitializer<Scalar>: Sized {
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self;

    fn new_fill(dimension: Index, fill: Scalar) -> Self
    where
        Scalar: Clone,
    {
        Self::new_func(dimension, |_, _| fill.clone())
    }

    fn from_matrix(matrix: &impl MatrixRef<Scalar>) -> Self
    where
        Scalar: Clone,
    {
        Self::new_func(matrix.dimension(), |row, column| {
            matrix.at(row, column).clone()
        })
    }
}

// Implementations for representation

impl<Scalar, Stored, Impl> MatrixRef<Scalar> for Representation<Stored, Impl>
where
    Stored: Borrow<Impl>,
    Impl: MatrixRef<Scalar>,
{
    fn dimension(&self) -> Index {
        self.represent().dimension()
    }

    fn at(&self, row: Index, column: Index) -> &Scalar {
        self.represent().at(row, column)
    }
}

impl<Scalar, Stored, Impl> MatrixMutRef<Scalar> for Representation<Stored, Impl>
where
    Stored: BorrowMut<Impl>,
    Impl: MatrixMutRef<Scalar>,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        self.represent_mut().at_mut(row, column)
    }
}

impl<Scalar, Stored, Impl> MatrixFuncInitializer<Scalar> for Representation<Stored, Impl>
where
    Stored: From<Impl>,
    Impl: MatrixFuncInitializer<Scalar>,
{
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self {
        Representation::from(Stored::from(Impl::new_func(dimension, fill)))
    }
}
