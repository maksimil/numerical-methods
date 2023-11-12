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

    pub fn new(dimension: Index, data: Vec<Scalar>) -> Self {
        Self { dimension, data }
    }
}

impl<Scalar> MatrixRef<Scalar> for DenseColMatrix<Scalar> {
    fn dimension(&self) -> Index {
        self.dimension
    }

    fn at(&self, row: Index, column: Index) -> &Scalar {
        &self.data[self.data_index(row, column)]
    }
}

impl<Scalar> MatrixMutRef<Scalar> for DenseColMatrix<Scalar> {
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        let data_index = self.data_index(row, column);
        &mut self.data[data_index]
    }
}

impl<Scalar> MatrixFuncInitializer<Scalar> for DenseColMatrix<Scalar> {
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self {
        let mut data = Vec::with_capacity(dimension * dimension);

        for column in 0..dimension {
            for row in 0..dimension {
                data.push(fill(row, column));
            }
        }

        DenseColMatrix::new(dimension, data)
    }
}

pub type DenseRowMatrix<Scalar> = super::transpose::MatrixTranspose<DenseColMatrix<Scalar>>;
