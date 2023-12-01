use std::{env::args, process::exit};

use iterative_methods::simple_iterative_solve;
use matrix::dense::DenseRowMatrix;
use qr_decomposition::QRDecomposition;
use test::Scalar;

use crate::{
    basic::Index,
    iterative_methods::zeidel_iterative_solve,
    lu_decomposition::LUDecomposition,
    matrix::{
        column::{ColumnFunc, ColumnFuncInitializer, ColumnRef},
        norms::NormedColumn,
    },
    test::{
        create_fifth_case, create_fifth_cases, create_static_test_cases, on_case,
        IterativeTestResult, TestResult,
    },
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

fn static_test_iterative_methods() {
    println!("Тест;bar x;e;МПИ;;;Метод Зейделя;;");
    println!(";;;x;d;k;x;d;k");

    for case in create_static_test_cases() {
        let test_n = case.name.clone();
        let precise_answer = case.answer.clone();
        println!("{test_n};{precise_answer:?};;;;;;;");
        for e in vec![1e-2, 1e-4, 1e-6, 1e-10] {
            let (it_answer, it_steps): (Vec<Scalar>, Index) =
                simple_iterative_solve(&case.matrix, &case.vector, e);
            let it_norm = ColumnFunc::new(precise_answer.dimension(), |i| {
                it_answer[i] - precise_answer[i]
            })
            .norm_one();

            let (zei_answer, zei_steps): (Vec<Scalar>, Index) =
                zeidel_iterative_solve(&case.matrix, &case.vector, e);
            let zei_norm = ColumnFunc::new(precise_answer.dimension(), |i| {
                zei_answer[i] - precise_answer[i]
            })
            .norm_one();

            println!(";;{e:.e};{it_answer:?};{it_norm:.e};{it_steps};{zei_answer:?};{zei_norm:.e};{zei_steps}");
        }
    }
}

fn dynamic_test_iterative_methods() {
    println!("Тест;n;epsilon;bar x;e;МПИ;;;Метод Зейделя;;");
    println!(";;;;;x;d;k;x;d;k");

    for case in create_fifth_cases(vec![1e-3, 1e-6], vec![4, 5, 8]) {
        let test_n = case.name.clone();
        let precise_answer = case.answer.clone();
        println!("{test_n};{precise_answer:?};;;;;;;");
        for e in vec![1e-2, 1e-4, 1e-6, 1e-10] {
            let (it_answer, it_steps): (Vec<Scalar>, Index) =
                simple_iterative_solve(&case.matrix, &case.vector, e);
            let it_norm = ColumnFunc::new(precise_answer.dimension(), |i| {
                it_answer[i] - precise_answer[i]
            })
            .norm_one();

            let (zei_answer, zei_steps): (Vec<Scalar>, Index) =
                zeidel_iterative_solve(&case.matrix, &case.vector, e);
            let zei_norm = ColumnFunc::new(precise_answer.dimension(), |i| {
                zei_answer[i] - precise_answer[i]
            })
            .norm_one();

            println!(";;;;{e:.e};{it_answer:?};{it_norm:.e};{it_steps};{zei_answer:?};{zei_norm:.e};{zei_steps}");
        }
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
        "static-iterative" => static_test_iterative_methods(),
        "dynamic-iterative" => dynamic_test_iterative_methods(),
        _ => {
            println!("Invalid test case");
            exit(1);
        }
    }
}

fn _main() {
    let matrix = DenseRowMatrix::new(3, vec![1.0, 0.0, 0.5, 0.9, 1.0, 0.0, 0.0, 0.0, 1.0]);
    let vector = vec![1.0, 0.0, 0.0];
    let (answer, iterations): (Vec<Scalar>, Index) = zeidel_iterative_solve(&matrix, &vector, 0.1);
    println!("Ans: {:?} steps={}", answer, iterations);
}
