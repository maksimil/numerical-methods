use std::borrow::{Borrow, BorrowMut};

use crate::basic::Index;

use super::traits::{MatrixFuncInitializer, MatrixMutRef, MatrixRef};

// Dense matrix in column major format
#[derive(Clone, Debug)]
pub struct DenseColMatrix<Scalar> {
    data: Vec<Scalar>,
    dimension: Index,
}

impl<Scalar> DenseColMatrix<Scalar> {
    fn data_index(&self, row: Index, column: Index) -> Index {
        row + column * self.dimension
    }
}

impl<Scalar, B> MatrixRef<Scalar> for B
where
    B: Borrow<DenseColMatrix<Scalar>>,
{
    fn dimension(&self) -> Index {
        self.borrow().dimension
    }

    fn at(&self, row: Index, column: Index) -> &Scalar {
        let data_index = self.borrow().data_index(row, column);
        &self.borrow().data[data_index]
    }
}

impl<Scalar, B> MatrixMutRef<Scalar> for B
where
    B: BorrowMut<DenseColMatrix<Scalar>>,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        let data_index = self.borrow_mut().data_index(row, column);
        &mut self.borrow_mut().data[data_index]
    }
}

impl<Scalar> MatrixFuncInitializer<Scalar> for DenseColMatrix<Scalar> {
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self {
        let mut data = Vec::with_capacity(dimension * dimension);

        for row in 0..dimension {
            for column in 0..dimension {
                data.push(fill(row, column));
            }
        }

        DenseColMatrix { data, dimension }
    }
}

pub type DenseRowMatrix<Scalar> = super::transpose::MatrixTranspose<DenseColMatrix<Scalar>>;
