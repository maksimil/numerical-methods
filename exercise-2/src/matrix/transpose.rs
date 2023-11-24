use crate::{
    basic::{Index, Numerical},
    representation::repr_ref,
};

use super::{
    column::{dot, ColumnOf},
    traits::{MatrixFuncInitializer, MatrixMutRef, MatrixRef},
};

// Represents a lazy transposed accessor for Matrix
#[derive(Clone, Debug)]
pub struct MatrixTranspose<Matrix> {
    matrix: Matrix,
}

impl<Matrix> From<Matrix> for MatrixTranspose<Matrix> {
    fn from(matrix: Matrix) -> Self {
        MatrixTranspose { matrix }
    }
}

impl<Matrix> MatrixRef for MatrixTranspose<Matrix>
where
    Matrix: MatrixRef,
{
    type Scalar = Matrix::Scalar;

    fn dimension(&self) -> Index {
        self.matrix.dimension()
    }

    fn at(&self, row: Index, column: Index) -> Matrix::Scalar {
        self.matrix.at(column, row)
    }
}

impl<Matrix> MatrixMutRef for MatrixTranspose<Matrix>
where
    Matrix: MatrixMutRef,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Matrix::Scalar {
        self.matrix.at_mut(column, row)
    }
}

impl<Matrix> MatrixFuncInitializer for MatrixTranspose<Matrix>
where
    Matrix: MatrixFuncInitializer,
{
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Matrix::Scalar) -> Self {
        Self::from(Matrix::new_func(dimension, |row, column| fill(column, row)))
    }
}

pub trait MatrixSymmetricSquare: MatrixRef {
    fn symmetric_square<OutMatrix>(&self) -> OutMatrix
    where
        OutMatrix: MatrixFuncInitializer + MatrixRef<Scalar = Self::Scalar>;
}

impl<T> MatrixSymmetricSquare for T
where
    T: MatrixRef,
    T::Scalar: Numerical,
{
    fn symmetric_square<OutMatrix>(&self) -> OutMatrix
    where
        OutMatrix: MatrixFuncInitializer + MatrixRef<Scalar = Self::Scalar>,
    {
        OutMatrix::new_func(self.dimension(), |i, j| {
            dot(
                &ColumnOf::new(repr_ref::<T>(&self), i),
                &ColumnOf::new(repr_ref::<T>(&self), j),
            )
        })
    }
}
