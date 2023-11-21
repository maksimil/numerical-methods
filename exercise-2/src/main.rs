use matrix::dense::DenseRowMatrix;

use crate::lu_decomposition::lu_decomposition::LUDecomposition;

mod basic;
mod lu_decomposition;
mod matrix;
mod representation;

type S = f64;

fn main() {
    let dimension = 4;
    let matrix = DenseRowMatrix::<S>::new(
        dimension,
        vec![
            2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
        ],
    );
    println!("A={:?}", matrix);

    let decomposition = LUDecomposition::calculate(matrix);
    println!("Decomposition={:#?}", decomposition);

    for index in 0..dimension {
        let mut v = vec![0.0; dimension];
        v[index] = 1.0;
        decomposition.solve(&mut v);
        println!("{}: {:?}", index, v);
    }
}
