use std::mem::swap;

use crate::basic::Numerical;

use super::traits::MatrixRef;

pub fn solve_upper<Matrix, Scalar>(matrix: &Matrix, vector: &mut Vec<Scalar>)
where
    Matrix: MatrixRef<Scalar>,
    Scalar: Numerical,
{
    let dimension = matrix.dimension();
    let mut b = vec![Scalar::zero(); dimension];
    swap(&mut b, vector);

    // // column major
    // for i in (0..dimension).rev() {
    //     vector[i] = vector[i] / matrix.at(i, i);
    //     for j in (0..i).rev() {
    //         vector[j] = vector[j] - matrix.at(j, i) * vector[i];
    //     }
    // }

    // row major
    for i in (0..dimension).rev() {
        vector[i] = (b[i]
            - ((i + 1)..dimension)
                .map(|j| matrix.at(i, j) * vector[j])
                .sum())
            / matrix.at(i, i);
    }
}
