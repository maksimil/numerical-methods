use std::borrow::{Borrow, BorrowMut};

use crate::{basic::Index, representation::Representation};

use super::traits::{MatrixMutRef, MatrixRef};

pub trait ColumnRef {
    type Scalar;
    fn dimension(&self) -> Index;
    fn at(&self, index: Index) -> Self::Scalar;
}

pub trait ColumnMut: ColumnRef {
    fn at_mut(&mut self, index: Index) -> &mut Self::Scalar;
}

pub trait ColumnInit: ColumnRef {
    fn new(dimension: Index, fill: Self::Scalar) -> Self;
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

impl<Scalar> ColumnInit for Vec<Scalar>
where
    Scalar: Clone,
{
    fn new(dimension: Index, fill: Self::Scalar) -> Self {
        vec![fill; dimension]
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
    pub fn of(matrix: Matrix, column: Index) -> Self {
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
