use rand_matrix::generate_matrix;

use crate::qr_algorithm::qr_algorithm;

mod inverse_power_method;
mod lu;
mod power_method;
mod qr_algorithm;
mod rand_matrix;
mod scalar;

fn main() {
    let dimension = 10;
    let generated = generate_matrix(dimension);

    println!("{:#?}", generated.eigenvalues);

    let matrix = generated.matrix;
    let res = qr_algorithm(&matrix, dimension);
    println!("res = {res:#?}");
}
