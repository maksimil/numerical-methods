use crate::basic::Index;

use super::traits::{MatrixMutRef, MatrixRef};

// Represents a lazy row permuted accessor for Matrix
#[derive(Clone, Debug)]
pub struct RowPermutedMatrix<Matrix> {
    matrix: Matrix,
    row_permutation: Vec<Index>,
}

impl<Matrix> RowPermutedMatrix<Matrix> {
    fn new(matrix: Matrix, row_permutation: Vec<Index>) -> RowPermutedMatrix<Matrix> {
        RowPermutedMatrix {
            matrix,
            row_permutation,
        }
    }
}

impl<Scalar, Matrix> MatrixRef<Scalar> for RowPermutedMatrix<Matrix>
where
    Matrix: MatrixRef<Scalar>,
{
    fn dimension(&self) -> Index {
        self.matrix.dimension()
    }

    fn at(&self, row: Index, column: Index) -> &Scalar {
        self.matrix.at(self.row_permutation[row], column)
    }
}

impl<Scalar, Matrix> MatrixMutRef<Scalar> for RowPermutedMatrix<Matrix>
where
    Matrix: MatrixMutRef<Scalar>,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        self.matrix.at_mut(self.row_permutation[row], column)
    }
}
