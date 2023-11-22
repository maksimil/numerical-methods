use std::mem::swap;

use crate::{
    basic::{Numerical, OtherNumericalOps},
    matrix::{column::ColumnMut, solve_upper::solve_upper, traits::MatrixMutRef},
};

fn sqrt<S>(value: S) -> S
where
    S: Numerical,
{
    value.sqrt_trait()
}

#[derive(Debug, Clone)]
pub struct QRDecomposition<Matrix>
where
    Matrix: MatrixMutRef,
    Matrix::Scalar: Numerical,
{
    // QA
    upper: Matrix,
    // hausdorf vectors in column major format, len=dim + (dim-1) + ...
    hausdorf_vectors: Vec<Matrix::Scalar>,
}

impl<Matrix> QRDecomposition<Matrix>
where
    Matrix: MatrixMutRef,
    Matrix::Scalar: Numerical,
{
    pub fn calculate(mut matrix: Matrix) -> QRDecomposition<Matrix> {
        let dimension = matrix.dimension();
        let mut hausdorf_vectors = vec![Matrix::Scalar::zero(); dimension * (dimension - 1)];

        for column in 0..dimension - 1 {
            let length = sqrt(
                (column..dimension)
                    .map(|row| matrix.at(row, column) * matrix.at(row, column))
                    .sum(),
            );

            // update the column and fill the hausdorf vector
            for row in column..dimension {
                swap(
                    &mut hausdorf_vectors[column * dimension + row],
                    matrix.at_mut(row, column),
                );
            }
            *matrix.at_mut(column, column) = length;
            hausdorf_vectors[column * dimension + column] -= length;

            let hausdorf_length = sqrt(
                (column..dimension)
                    .map(|row| {
                        hausdorf_vectors[column * dimension + row]
                            * hausdorf_vectors[column * dimension + row]
                    })
                    .sum(),
            );

            for row in column..dimension {
                hausdorf_vectors[column * dimension + row] /= hausdorf_length;
            }

            // updating the rest of the matrix
            for affected_column in column + 1..dimension {
                let dot_product = (column..dimension)
                    .map(|i| {
                        hausdorf_vectors[column * dimension + i] * matrix.at(i, affected_column)
                    })
                    .sum();
                for row in column..dimension {
                    *matrix.at_mut(row, affected_column) -=
                        (hausdorf_vectors[column * dimension + row] * dot_product)
                            * Matrix::Scalar::from(2);
                }
            }
        }

        Self {
            upper: matrix,
            hausdorf_vectors,
        }
    }

    pub fn solve<Column>(&self, vector: &mut Column)
    where
        Column: ColumnMut<Scalar = Matrix::Scalar>,
    {
        let dimension = self.upper.dimension();

        // applying Q to vector
        for hausdorf_index in 0..dimension - 1 {
            let dot_product = (hausdorf_index..dimension)
                .map(|i| self.hausdorf_vectors[hausdorf_index * dimension + i] * vector.at(i))
                .sum();
            for row in hausdorf_index..dimension {
                *vector.at_mut(row) -= (self.hausdorf_vectors[hausdorf_index * dimension + row]
                    * dot_product)
                    * Matrix::Scalar::from(2);
            }
        }

        // solving upper system
        solve_upper(&self.upper, vector);
    }
}
