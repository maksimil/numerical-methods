use crate::basic::Index;

use super::{
    permutation::Permutation,
    traits::{MatrixMutRef, MatrixRef},
};

// Represents a lazy row permuted accessor for P^T M
#[derive(Clone, Debug)]
pub struct RowPermutedMatrix<Matrix, Perm> {
    matrix: Matrix,
    row_permutation: Perm,
}

impl<Matrix, Perm> RowPermutedMatrix<Matrix, Perm> {
    pub fn new(matrix: Matrix, row_permutation: Perm) -> RowPermutedMatrix<Matrix, Perm> {
        RowPermutedMatrix {
            matrix,
            row_permutation,
        }
    }
}

impl<Matrix, Perm> MatrixRef for RowPermutedMatrix<Matrix, Perm>
where
    Matrix: MatrixRef,
    Perm: Permutation,
{
    type Scalar = Matrix::Scalar;

    fn dimension(&self) -> Index {
        self.matrix.dimension()
    }

    fn at(&self, row: Index, column: Index) -> Self::Scalar {
        self.matrix.at(self.row_permutation.permute(row), column)
    }
}

impl<Matrix, Perm> MatrixMutRef for RowPermutedMatrix<Matrix, Perm>
where
    Matrix: MatrixMutRef,
    Perm: Permutation,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Self::Scalar {
        self.matrix
            .at_mut(self.row_permutation.permute(row), column)
    }
}
