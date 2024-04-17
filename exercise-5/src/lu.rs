use std::mem::swap;

use crate::scalar::*;

pub fn dot_product(x: &[Scalar], y: &[Scalar]) -> Scalar {
    debug_assert!(x.len() == y.len());

    let mut r = SCALAR_ZERO;

    for i in 0..x.len() {
        r += x[i] * y[i];
    }

    r
}

const UPPER: Index = 0;
const LOWER: Index = 1;

fn triangular_solve<const TRIANGULARITY: Index>(
    matrix: &Vec<Scalar>,
    permutation: &Vec<Index>,
    x: &mut Vec<Scalar>,
) {
    assert!(TRIANGULARITY == UPPER || TRIANGULARITY == LOWER);
    let dimension = permutation.len();

    for kk in 0..dimension {
        let k = if TRIANGULARITY == UPPER {
            permutation[dimension - 1 - kk]
        } else {
            permutation[kk]
        };

        let v = x[k];
        x[k] = SCALAR_ZERO;
        x[k] = (v - dot_product(&matrix[k * dimension..k * dimension + dimension], &x))
            / matrix[k * dimension + k];
    }
}

#[derive(Debug, Clone)]
pub struct LUDecomposition {
    // matricies are in row-major format
    lower_: Vec<Scalar>,
    upper_: Vec<Scalar>,

    permutation_: Vec<Index>,
    inverse_permutation_: Vec<Index>,
}

const MIN_PIVOT: Scalar = 1e-7;

impl LUDecomposition {
    pub fn new_zero(dimension: Index) -> LUDecomposition {
        let lower = vec![SCALAR_ZERO; dimension * dimension];
        let upper = vec![SCALAR_ZERO; dimension * dimension];

        let permutation = vec![0; dimension];
        let inverse_permutation = vec![0; dimension];

        LUDecomposition {
            lower_: lower,
            upper_: upper,
            permutation_: permutation,
            inverse_permutation_: inverse_permutation,
        }
    }

    // input matrix is in row-major format
    pub fn compute(matrix: &Vec<Scalar>, dimension: Index) -> Option<LUDecomposition> {
        let mut decomposition = LUDecomposition::new_zero(dimension);
        let mut row_used = vec![false; dimension];

        // setting permutations to identity
        for k in 0..dimension {
            decomposition.permutation_[k] = k;
            decomposition.inverse_permutation_[k] = k;
        }

        // setting L to identity
        for k in 0..dimension {
            decomposition.lower_[dimension * k + k] = SCALAR_ONE;
        }

        // computing LU
        for k in 0..dimension {
            // computing the column
            let mut column = vec![SCALAR_ZERO; dimension];

            for i in 0..dimension {
                column[i] = matrix[dimension * i + k];
            }

            triangular_solve::<LOWER>(
                &decomposition.lower_,
                &decomposition.permutation_,
                &mut column,
            );

            // choosing pivot
            let pivot_row = {
                let mut pivot_row = 0;

                for i in 0..dimension {
                    if !row_used[i]
                        && (row_used[pivot_row] || column[i].abs() > column[pivot_row].abs())
                    {
                        pivot_row = i;
                    }
                }

                pivot_row
            };

            let pivot_value = column[pivot_row];

            if pivot_value.abs() < MIN_PIVOT {
                return None;
            }

            row_used[pivot_row] = true;

            // computing L and U
            for i in 0..dimension {
                if row_used[i] {
                    decomposition.upper_[dimension * i + pivot_row] = column[i];
                } else {
                    decomposition.lower_[dimension * i + pivot_row] = column[i] / pivot_value;
                }
            }

            // computing the permutation
            let ipr = decomposition.inverse_permutation_[pivot_row];
            let pk = decomposition.permutation_[k];

            decomposition.permutation_[k] = pivot_row;
            decomposition.inverse_permutation_[pivot_row] = k;
            decomposition.permutation_[ipr] = pk;
            decomposition.inverse_permutation_[pk] = ipr;
        }

        Some(decomposition)
    }

    pub fn solve(&self, x: &mut Vec<Scalar>) {
        triangular_solve::<LOWER>(&self.lower_, &self.permutation_, x);
        triangular_solve::<UPPER>(&self.upper_, &self.permutation_, x);

        let mut temp = vec![SCALAR_ZERO; x.len()];

        for i in 0..x.len() {
            temp[i] = x[self.permutation_[i]];
        }

        swap(&mut temp, x);
    }
}
