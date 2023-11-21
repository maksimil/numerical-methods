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

impl<Scalar, Matrix, Perm> MatrixRef<Scalar> for RowPermutedMatrix<Matrix, Perm>
where
    Matrix: MatrixRef<Scalar>,
    Perm: Permutation,
{
    fn dimension(&self) -> Index {
        self.matrix.dimension()
    }

    fn at(&self, row: Index, column: Index) -> Scalar {
        self.matrix.at(self.row_permutation.permute(row), column)
    }
}

impl<Scalar, Matrix, Perm> MatrixMutRef<Scalar> for RowPermutedMatrix<Matrix, Perm>
where
    Matrix: MatrixMutRef<Scalar>,
    Perm: Permutation,
{
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar {
        self.matrix
            .at_mut(self.row_permutation.permute(row), column)
    }
}
