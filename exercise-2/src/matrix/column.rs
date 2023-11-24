use std::borrow::{Borrow, BorrowMut};

use crate::{
    basic::{Index, Numerical, OtherNumericalOps},
    representation::{repr_ref, Representation},
};

use super::{
    traits::{MatrixMutRef, MatrixRef},
    transpose::MatrixTranspose,
};

pub trait ColumnRef {
    type Scalar;
    fn dimension(&self) -> Index;
    fn at(&self, index: Index) -> Self::Scalar;
}

pub trait ColumnMut: ColumnRef {
    fn at_mut(&mut self, index: Index) -> &mut Self::Scalar;
}

pub trait ColumnFuncInitializer: Sized + ColumnRef {
    fn new_func(dimension: Index, fill: impl Fn(Index) -> Self::Scalar) -> Self;

    fn new_fill(dimension: Index, fill: Self::Scalar) -> Self
    where
        Self::Scalar: Clone,
    {
        Self::new_func(dimension, |_| fill.clone())
    }

    fn from_column<S>(column: &impl ColumnRef<Scalar = S>) -> Self
    where
        S: Into<Self::Scalar>,
    {
        Self::new_func(column.dimension(), |i| column.at(i).into())
    }
}

impl<Scalar> ColumnRef for Vec<Scalar>
where
    Scalar: Clone,
{
    type Scalar = Scalar;

    fn dimension(&self) -> Index {
        self.len()
    }

    fn at(&self, index: Index) -> Scalar {
        self[index].clone()
    }
}

impl<Scalar> ColumnMut for Vec<Scalar>
where
    Scalar: Clone,
{
    fn at_mut(&mut self, index: Index) -> &mut Scalar {
        &mut self[index]
    }
}

impl<Scalar> ColumnFuncInitializer for Vec<Scalar>
where
    Scalar: Clone,
{
    fn new_func(dimension: Index, fill: impl Fn(Index) -> Self::Scalar) -> Self {
        let mut data = Vec::with_capacity(dimension);
        for i in 0..dimension {
            data.push(fill(i));
        }
        data
    }
}

pub fn dot<Rhs, Lhs>(rhs: &Rhs, lhs: &Lhs) -> Rhs::Scalar
where
    Rhs: ColumnRef,
    Lhs: ColumnRef<Scalar = Rhs::Scalar>,
    Rhs::Scalar: Numerical,
{
    (0..rhs.dimension()).map(|i| rhs.at(i) * lhs.at(i)).sum()
}

pub fn apply<Matrix, ColumnIn, ColumnOut>(matrix: &Matrix, column: &ColumnIn) -> ColumnOut
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnFuncInitializer + ColumnRef<Scalar = Matrix::Scalar>,
{
    ColumnOut::new_func(matrix.dimension(), |i| {
        dot(
            &ColumnOf::new(MatrixTranspose::from(repr_ref::<Matrix>(matrix)), i),
            column,
        )
    })
}

pub struct ColumnFunc<Scalar, F>
where
    F: Fn(Index) -> Scalar,
{
    dimension: Index,
    function: F,
}

impl<Scalar, F> ColumnFunc<Scalar, F>
where
    F: Fn(Index) -> Scalar,
{
    pub fn new(dimension: Index, function: F) -> Self {
        Self {
            dimension,
            function,
        }
    }
}

impl<Scalar, F> ColumnRef for ColumnFunc<Scalar, F>
where
    F: Fn(Index) -> Scalar,
{
    type Scalar = Scalar;

    fn dimension(&self) -> Index {
        self.dimension
    }

    fn at(&self, index: Index) -> Self::Scalar {
        (self.function)(index)
    }
}

pub struct ColumnOf<Matrix>
where
    Matrix: MatrixRef,
{
    matrix: Matrix,
    column: Index,
}

impl<Matrix> ColumnOf<Matrix>
where
    Matrix: MatrixRef,
{
    pub fn new(matrix: Matrix, column: Index) -> Self {
        Self { matrix, column }
    }
}

impl<Matrix> ColumnRef for ColumnOf<Matrix>
where
    Matrix: MatrixRef,
{
    type Scalar = Matrix::Scalar;

    fn dimension(&self) -> Index {
        self.matrix.dimension()
    }

    fn at(&self, index: Index) -> Matrix::Scalar {
        self.matrix.at(index, self.column)
    }
}

impl<Matrix> ColumnMut for ColumnOf<Matrix>
where
    Matrix: MatrixMutRef,
{
    fn at_mut(&mut self, index: Index) -> &mut Matrix::Scalar {
        self.matrix.at_mut(index, self.column)
    }
}

impl<Stored, Impl> ColumnRef for Representation<Stored, Impl>
where
    Stored: Borrow<Impl>,
    Impl: ColumnRef,
{
    type Scalar = Impl::Scalar;

    fn dimension(&self) -> Index {
        self.represent().dimension()
    }

    fn at(&self, index: Index) -> Self::Scalar {
        self.represent().at(index)
    }
}

impl<Stored, Impl> ColumnMut for Representation<Stored, Impl>
where
    Stored: BorrowMut<Impl>,
    Impl: ColumnMut,
{
    fn at_mut(&mut self, index: Index) -> &mut Self::Scalar {
        self.represent_mut().at_mut(index)
    }
}
