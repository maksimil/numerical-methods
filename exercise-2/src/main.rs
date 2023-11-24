use std::{env::args, process::exit};

use matrix::dense::DenseRowMatrix;
use qr_decomposition::QRDecomposition;
use test::Scalar;

use crate::{
    lu_decomposition::LUDecomposition,
    test::{create_fifth_case, create_fifth_cases, create_static_test_cases, on_case, TestResult},
};

mod basic;
mod iterative_methods;
mod lu_decomposition;
mod matrix;
mod qr_decomposition;
mod representation;
mod test;

fn solve_lu(matrix: &DenseRowMatrix<Scalar>, vector: &Vec<Scalar>) -> Vec<Scalar> {
    let decomposition = LUDecomposition::calculate(matrix.clone());
    let mut vector_mut = vector.clone();
    decomposition.solve(&mut vector_mut);
    vector_mut
}

fn solve_qr(matrix: &DenseRowMatrix<Scalar>, vector: &Vec<Scalar>) -> Vec<Scalar> {
    let decomposition = QRDecomposition::calculate(matrix.clone());
    let mut vector_mut = vector.clone();
    decomposition.solve(&mut vector_mut);
    vector_mut
}

fn static_test_direct_methods() {
    println!("Тест;bar x;LU;;QR;");
    println!(";;x;d;x;d");

    for case in create_static_test_cases() {
        let test_n = case.name.clone();
        let precise_answer = case.answer.clone();
        let TestResult {
            answer: lu_answer,
            norm_one: lu_error,
        } = on_case(&case, solve_lu);
        let TestResult {
            answer: qr_answer,
            norm_one: qr_error,
        } = on_case(&case, solve_qr);

        println!("{test_n};{precise_answer:?};{lu_answer:?};{lu_error};{qr_answer:?};{qr_error}");
    }
}

fn dynamic_test_direct_methods() {
    println!("Тест;n;epsilon;bar x;LU;;QR;");
    println!(";;;;x;d;x;d");

    for case in create_fifth_cases(vec![1e-3, 1e-6, 1e-9, 1e-12], vec![4, 5, 8, 16, 32]) {
        let test_n = case.name.clone();
        let precise_answer = case.answer.clone();
        let TestResult {
            answer: lu_answer,
            norm_one: lu_error,
        } = on_case(&case, solve_lu);
        let TestResult {
            answer: qr_answer,
            norm_one: qr_error,
        } = on_case(&case, solve_qr);

        println!("{test_n};{precise_answer:?};{lu_answer:?};{lu_error};{qr_answer:?};{qr_error}");
    }
}

fn main() {
    let cli_args: Vec<String> = args().collect();

    if cli_args.len() != 2 {
        println!(
            "usage: [command] [static-direct|dynamic-direct|static-iterative|dynamic-iterative]"
        );
        exit(1);
    }

    let test_type = cli_args[1].as_str();

    match test_type {
        "static-direct" => static_test_direct_methods(),
        "dynamic-direct" => dynamic_test_direct_methods(),
        _ => {
            println!("Invalid test case");
            exit(1);
        }
    }
}
