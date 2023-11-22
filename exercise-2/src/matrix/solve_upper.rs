use crate::basic::Numerical;

use super::{column::ColumnMut, traits::MatrixRef};

pub fn solve_upper<Matrix, Column>(matrix: &Matrix, vector: &mut Column)
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    Column: ColumnMut<Scalar = Matrix::Scalar>,
{
    let dimension = matrix.dimension();

    // column major
    for i in (0..dimension).rev() {
        *vector.at_mut(i) = vector.at(i) / matrix.at(i, i);
        for j in (0..i).rev() {
            *vector.at_mut(j) = vector.at(j) - matrix.at(j, i) * vector.at(i);
        }
    }

    // // row major
    // for i in (0..dimension).rev() {
    //     vector[i] = (b[i]
    //         - ((i + 1)..dimension)
    //             .map(|j| matrix.at(i, j) * vector[j])
    //             .sum())
    //         / matrix.at(i, i);
    // }
}
