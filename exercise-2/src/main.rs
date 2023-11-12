use matrix::dense::DenseColMatrix;

use crate::{
    matrix::{
        traits::{MatrixFuncInitializer, MatrixRef},
        transpose::MatrixTranspose,
    },
    representation::repr,
};

mod basic;
mod matrix;
mod representation;

fn main() {
    let col_matrix = DenseColMatrix::new(3, vec![1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 3.0, 2.0, 1.0]);
    let transpose_lazy = MatrixTranspose::from(repr(&col_matrix));
    let transpose_calculated = DenseColMatrix::from_matrix(&transpose_lazy);
    println!("A[0,1]={:?}", col_matrix.at(0, 1));
    println!("A={:?}", col_matrix);
    println!("lazy A^T={:?}", transpose_lazy);
    println!("A^T={:?}", transpose_calculated);
}
