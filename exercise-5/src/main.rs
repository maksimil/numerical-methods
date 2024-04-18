use rand::distributions::Distribution;
use rand_matrix::generate_matrix;

use crate::inverse_power_method::inverse_power_method;
use crate::power_method::power_method;
use crate::qr_algorithm::qr_algorithm;
use crate::rand_matrix::{rand_scalar, GenerateMatrix};
use crate::scalar::*;

mod inverse_power_method;
mod lu;
mod power_method;
mod qr_algorithm;
mod rand_matrix;
mod scalar;

fn main() {
    // read dimension
    let dimension = {
        let mut buf = String::new();
        let res = std::io::stdin().read_line(&mut buf);

        if let Err(e) = res {
            eprintln!("Failed to parse dimension: {e}");
        }

        let r = buf.trim().parse::<Index>();

        match r {
            Ok(dimension) => dimension,
            Err(e) => {
                eprintln!("Failed to parse dimension: {e}");
                std::process::exit(1);
            }
        }
    };

    // generate matrix
    let GenerateMatrix {
        matrix,
        eigenvalues,
    } = generate_matrix(dimension);

    // power method
    println!("Power method");
    let pm_result = power_method(&matrix, dimension);
    let expected_result = if eigenvalues[dimension - 1].abs() > eigenvalues[0].abs() {
        eigenvalues[dimension - 1]
    } else {
        eigenvalues[0]
    };

    match pm_result {
        Some(pm) => {
            let eigenvalue = pm.value;
            let error = (eigenvalue - expected_result).abs();
            println!("Expected = {expected_result}");
            println!("Result = {eigenvalue}");
            println!("Error = {error:.e}");
        }
        None => {
            println!("Did not converge");
        }
    }
    println!();

    // inverse power method
    println!("Inverse power method");
    let eigen_idx =
        rand::distributions::Uniform::from(0..dimension).sample(&mut rand::thread_rng());
    let expected_result = eigenvalues[eigen_idx];
    let approximation = expected_result + rand_scalar(VALUES_SEPARATION * 0.5);
    let ipm_result = inverse_power_method(&matrix, dimension, approximation);

    match ipm_result {
        Some(ipm) => {
            let eigenvalue = ipm.value;
            let error = (eigenvalue - expected_result).abs();
            let init_error = (eigenvalue - approximation).abs();
            println!("Expected = {expected_result}");
            println!("Initial estimate = {approximation}");
            println!("Initial error = {init_error:.e}");
            println!("Result = {eigenvalue}");
            println!("Error = {error:.e}");
        }
        None => {
            println!("Did not converge");
        }
    }

    // qr algorithm
    println!("QR-algorithm");
    let qr_result = qr_algorithm(&matrix, dimension);

    match qr_result {
        Some(values) => {
            let mut average_error = SCALAR_ZERO;

            for k in 0..dimension {
                let expected = eigenvalues[k];
                let calculated = values[k];
                let error = (expected - calculated).abs();
                average_error += error;
                println!(
                    "Error = {:>9} Expected = {:>9} Result = {:>9}",
                    format!("{error:.3e}"),
                    format!("{expected:.3}"),
                    format!("{calculated:.3}")
                );
            }
            average_error /= dimension as Scalar;
            println!("Average error = {average_error:.5e}");
        }
        None => {
            println!("Did not converge");
        }
    }
}
