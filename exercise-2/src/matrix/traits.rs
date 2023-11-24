use std::borrow::{Borrow, BorrowMut};

use crate::{basic::Index, representation::Representation};

pub trait MatrixRef {
    type Scalar;

    fn dimension(&self) -> Index;
    fn at(&self, row: Index, column: Index) -> Self::Scalar;
}

pub trait MatrixMutRef: MatrixRef {
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Self::Scalar;
}

pub trait MatrixFuncInitializer: Sized + MatrixRef {
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Self::Scalar) -> Self;

    fn new_fill(dimension: Index, fill: Self::Scalar) -> Self
    where
        Self::Scalar: Clone,
    {
        Self::new_func(dimension, |_, _| fill.clone())
    }

    fn from_matrix<S>(matrix: &impl MatrixRef<Scalar = S>) -> Self
    where
        S: Into<Self::Scalar>,
    {
        Self::new_func(matrix.dimension(), |row, column| {
            matrix.at(row, column).into()
        })
    }
}

pub struct MatrixFunc<Scalar, F: Fn(Index, Index) -> Scalar> {
    dimension: Index,
    function: F,
}

impl<Scalar, F> MatrixFunc<Scalar, F>
where
    F: Fn(Index, Index) -> Scalar,
{
    pub fn new(dimension: Index, function: F) -> Self {
        Self {
            dimension,
            function,
        }
    }
}

impl<Scalar, F> MatrixRef for MatrixFunc<Scalar, F>
where
    F: Fn(Index, Index) -> Scalar,
{
    type Scalar = Scalar;

    fn dimension(&self) -> Index {
        self.dimension
    }

    fn at(&self, row: Index, column: Index) -> Self::Scalar {
        (self.function)(row, column)
    }
}

// Implementations for representation

impl<Stored, Impl> MatrixRef for Representation<Stored, Impl>
where
    Stored: Borrow<Impl>,
    Impl: MatrixRef,
    Impl::Scalar: Clone,
{
    type Scalar = Impl::Scalar;

    fn dimension(&self) -> Index {
        self.represent().dimension()
    }

    fn at(&self, row: Index, column: Index) -> Impl::Scalar {
        self.represent().at(row, column)
    }
}

impl<Stored, Impl> MatrixMutRef for Representation<Stored, Impl>
where
    Stored: BorrowMut<Impl>,
    Impl: MatrixMutRef,
    Impl::Scalar: Clone,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Impl::Scalar {
        self.represent_mut().at_mut(row, column)
    }
}

impl<Stored, Impl> MatrixFuncInitializer for Representation<Stored, Impl>
where
    Stored: Borrow<Impl> + From<Impl>,
    Impl: MatrixFuncInitializer,
    Impl::Scalar: Clone,
{
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Impl::Scalar) -> Self {
        Representation::from(Stored::from(Impl::new_func(dimension, fill)))
    }
}
