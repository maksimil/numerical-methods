use std::mem::swap;

use crate::{
    basic::{Index, Numerical, INDEX_NOT_FOUND},
    matrix::{
        row_permuted::RowPermutedMatrix,
        solve_upper::solve_upper,
        traits::{MatrixMutRef, MatrixRef},
    },
    representation::repr_ref,
};

#[derive(Debug, Clone)]
pub struct LUDecomposition<Matrix, Scalar> {
    // L'A
    upper_permuted: Matrix,
    // P such that PL'A is upper
    permutation: Vec<Index>,
    inverse_permutation: Vec<Index>,
    // eta vectors in row major format, len=dim*(dim-1),
    eta_data: Vec<Scalar>,
}

impl<Matrix, Scalar> LUDecomposition<Matrix, Scalar> {
    pub fn calculate(mut matrix: Matrix) -> LUDecomposition<Matrix, Scalar>
    where
        Scalar: Numerical,
        Matrix: MatrixRef<Scalar> + MatrixMutRef<Scalar>,
    {
        let dimension = matrix.dimension();
        let mut permutation = vec![INDEX_NOT_FOUND; dimension];
        let mut inverse_permutation = vec![INDEX_NOT_FOUND; dimension];
        let mut eta_data = vec![Scalar::zero(); dimension * (dimension - 1)];

        for column in 0..dimension - 1 {
            let (_, pivot_row) = (0..dimension)
                .filter(|row| permutation[*row] == INDEX_NOT_FOUND)
                .fold(
                    (Scalar::zero(), INDEX_NOT_FOUND),
                    |(max_value, max_index), row| {
                        let row_value = matrix.at(row, column).abs_trait();
                        if row_value >= max_value {
                            (row_value, row)
                        } else {
                            (max_value, max_index)
                        }
                    },
                );

            println!("pivot={},{}", pivot_row, column);

            permutation[pivot_row] = column;
            inverse_permutation[column] = pivot_row;

            for row in 0..dimension {
                if permutation[row] == INDEX_NOT_FOUND {
                    let factor = -matrix.at(row, column) / matrix.at(pivot_row, column);
                    eta_data[column * dimension + row] = factor.clone();

                    *matrix.at_mut(row, column) = Scalar::zero();

                    if factor != Scalar::zero() {
                        for affected_column in column + 1..dimension {
                            let fill_in = matrix.at(row, affected_column)
                                + factor.clone() * matrix.at(pivot_row, affected_column);
                            *matrix.at_mut(row, affected_column) = fill_in;
                        }
                    }
                }
            }
        }

        for row in 0..dimension {
            if permutation[row] == INDEX_NOT_FOUND {
                permutation[row] = dimension - 1;
                inverse_permutation[dimension - 1] = row;
            }
        }

        LUDecomposition {
            permutation,
            inverse_permutation,
            eta_data,
            upper_permuted: matrix,
        }
    }

    pub fn solve(&self, vector: &mut Vec<Scalar>)
    where
        Scalar: Numerical,
        Matrix: MatrixRef<Scalar>,
    {
        let dimension = self.upper_permuted.dimension();
        // applying lower triangular matricies
        for eta_index in 0..dimension - 1 {
            let eta_column = self.inverse_permutation[eta_index];
            let vector_value = vector[eta_column];
            for index in 0..dimension {
                vector[index] =
                    vector[index] + self.eta_data[eta_index * dimension + index] * vector_value;
            }
        }

        // permutation
        {
            let mut unpermuted = vec![Scalar::zero(); dimension];
            swap(&mut unpermuted, vector);
            for index in 0..dimension {
                vector[index] = unpermuted[self.inverse_permutation[index]];
            }
        }

        // solving upper triangular system
        solve_upper(
            &RowPermutedMatrix::new(
                repr_ref(&self.upper_permuted),
                repr_ref(&self.inverse_permutation),
            ),
            vector,
        );
    }
}
