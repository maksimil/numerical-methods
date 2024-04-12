use rand_matrix::generate_matrix;

use crate::power_method::power_method;

mod lu;
mod power_method;
mod rand_matrix;
mod scalar;

fn main() {
    let dimension = 100;
    let generated = generate_matrix(dimension);

    println!("{:#?}", generated.eigenvalues);

    let matrix = generated.matrix;

    let pm = power_method(&matrix, dimension);

    println!("{pm:#?}");
}
