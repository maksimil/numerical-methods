use std::borrow::{Borrow, BorrowMut};

use crate::{basic::Index, representation::Representation};

use super::traits::{MatrixFuncInitializer, MatrixMutRef, MatrixRef};

#[derive(Clone, Debug)]
pub struct MatrixTranspose<MatrixStored, MatrixImpl> {
    matrix: Representation<MatrixStored, MatrixImpl>,
}

impl<MatrixStored, MatrixImpl> From<MatrixStored> for MatrixTranspose<MatrixStored, MatrixImpl> {
    pub fn from(matrix: MatrixStored) -> Self {
        MatrixTranspose {
            matrix: Representation::from(matrix),
        }
    }
}

impl<Scalar, MatrixStored, MatrixImpl> MatrixRef<Scalar>
    for MatrixTranspose<MatrixStored, MatrixImpl>
where
    MatrixImpl: MatrixRef<Scalar>,
    MatrixStored: Borrow<MatrixImpl>,
{
    fn dimension(&self) -> Index {
        self.matrix.represent().dimension()
    }

    fn at(&self, row: Index, column: Index) -> &Scalar {
        self.matrix.represent().at(column, row)
    }
}

impl<Scalar, MatrixStored, MatrixImpl> MatrixMutRef<Scalar>
    for MatrixTranspose<MatrixStored, MatrixImpl>
where
    MatrixImpl: MatrixMutRef<Scalar>,
    MatrixStored: BorrowMut<MatrixImpl>,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        self.matrix.represent_mut().at_mut(column, row)
    }
}

impl<Scalar, MatrixStored, MatrixImpl> MatrixFuncInitializer<Scalar>
    for MatrixTranspose<MatrixStored, MatrixImpl>
where
    MatrixStored: MatrixFuncInitializer<Scalar>,
{
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self {
        Self::from(MatrixStored::new_func(dimension, |row, column| {
            fill(column, row)
        }))
    }
}
