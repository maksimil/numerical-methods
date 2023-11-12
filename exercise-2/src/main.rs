use crate::matrix::{
    dense::DenseColMatrix,
    traits::{MatrixFillInitializer, MatrixMutRef},
};

mod basic;
mod matrix;

fn main() {
    let mut matrix = DenseColMatrix::new_fill(2, 0.0);
    *matrix.at_mut(0, 1) = 1.0;
    *matrix.at_mut(0, 0) = 1.0;
    *matrix.at_mut(1, 0) = -1.0;
    println!("A: {:?}", matrix);
    // println!("A^TA: {:?}", matrix.transpose_mul());
    // println!(
    //     "norm_inf: {}, norm_1: {}",
    //     matrix.norm_inf(),
    //     matrix.norm_1()
    // );
}
