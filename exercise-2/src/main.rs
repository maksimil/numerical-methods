use matrix::{
    dense::DenseColMatrix,
    traits::{MatrixFillInitializer, MatrixMutRef, MatrixRef},
    transpose::MatrixTranspose,
};

mod basic;
mod matrix;

fn main() {
    let mut matrix = DenseColMatrix::new_fill(2, 0.0 as f64);
    *matrix.at_mut(0, 1) = 1.0;
    *matrix.at_mut(0, 0) = 1.0;
    *matrix.at_mut(1, 0) = -1.0;

    println!("A: {:?}", matrix);

    let mut lazy_transpose = MatrixTranspose::from_matrix(&mut matrix);

    println!("A^T[0,1]={:?}", lazy_transpose.at(0, 1));

    *lazy_transpose.at_mut(0, 1) = 2.0;

    println!("A: {:?}", matrix);
    // println!("A^TA: {:?}", matrix.transpose_mul());
    // println!(
    //     "norm_inf: {}, norm_1: {}",
    //     matrix.norm_inf(),
    //     matrix.norm_1()
    // );
}
