use crate::basic::Index;

use super::traits::{MatrixFuncInitializer, MatrixMutRef, MatrixRef};

#[derive(Clone, Debug)]
pub struct MatrixTranspose<Matrix> {
    matrix: Matrix,
}

impl<Matrix> MatrixTranspose<Matrix> {
    pub fn from_matrix(matrix: Matrix) -> Self {
        MatrixTranspose { matrix }
    }
}

impl<Scalar, Matrix> MatrixRef<Scalar> for MatrixTranspose<Matrix>
where
    Matrix: MatrixRef<Scalar>,
{
    fn dimension(&self) -> Index {
        self.matrix.dimension()
    }

    fn at(&self, row: Index, column: Index) -> &Scalar {
        self.matrix.at(column, row)
    }
}

impl<Scalar, Matrix> MatrixMutRef<Scalar> for MatrixTranspose<Matrix>
where
    Matrix: MatrixMutRef<Scalar>,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        self.matrix.at_mut(column, row)
    }
}

impl<Scalar, Matrix> MatrixFuncInitializer<Scalar> for MatrixTranspose<Matrix>
where
    Matrix: MatrixFuncInitializer<Scalar>,
{
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self {
        Self::from_matrix(Matrix::new_func(dimension, |row, column| fill(column, row)))
    }
}
