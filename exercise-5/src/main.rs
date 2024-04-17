use rand::distributions::{Distribution, Uniform};
use rand_matrix::generate_matrix;

use crate::{
    inverse_power_method::inverse_power_method,
    rand_matrix::{rand_scalar, VALUES_SEPARATION},
};

mod inverse_power_method;
mod lu;
mod power_method;
mod rand_matrix;
mod scalar;

fn main() {
    let dimension = 1000;
    let generated = generate_matrix(dimension);

    println!("{:#?}", generated.eigenvalues);

    let matrix = generated.matrix;

    let value_index = Uniform::from(0..dimension).sample(&mut rand::thread_rng());
    let approximation = generated.eigenvalues[value_index] + rand_scalar(VALUES_SEPARATION / 2.);

    let pm = inverse_power_method(&matrix, dimension, approximation);

    println!("approximation = {approximation:?}");
    println!("pm = {pm:#?}");
}
