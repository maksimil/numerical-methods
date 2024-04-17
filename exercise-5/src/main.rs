use rand_matrix::generate_matrix;

use crate::qr_algorithm::hessenberg_form;

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
    let mut hessenberg = matrix.clone();

    hessenberg_form(&mut hessenberg, dimension);

    println!("A = {matrix:#?}");
    println!("B = {hessenberg:#?}");
}
