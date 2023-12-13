use basic::{Interval, Scalar};
use exercise_2::matrix::dense::DenseRowMatrix;
use newton_method::newton_method;
use newton_method_n::newton_method_n;

use crate::newton_method_n::newton_method_n_interpolated;

pub mod basic;
pub mod localize_root;
pub mod newton_method;
pub mod newton_method_n;

fn task_f(x: Scalar) -> Scalar {
    (Scalar::from(0.5) * x + Scalar::from(0.2)).tan() - x * x
}

fn task_fprime(x: Scalar) -> Scalar {
    Scalar::from(1.0)
        / (Scalar::from(2.0) * (x * Scalar::from(0.5) + Scalar::from(0.2)).cos().powi(2))
        - Scalar::from(2.0) * x
}

fn test_newton_method() {
    println!("\x1b[32m=== Scalar newton method ===\x1b[0m");
    let interval = Interval::new(0.0, 9.0);
    let max_iterations = 100;
    let answer = newton_method(task_f, task_fprime, interval, 1e-4, max_iterations);
    match answer {
        Some(x) => {
            println!("Newton method result = {x}");
            println!("Residual = {r:.e}", r = task_f(x));
        }
        None => println!("Newton method did not terminate in {max_iterations}"),
    }
}

fn task_fn(v: &Vec<Scalar>) -> Vec<Scalar> {
    let x = v[0];
    let y = v[1];
    vec![-x + (y + 0.5).sin() - 1.0, (x - 2.0).cos() + y]
}

fn task_fnprime(v: &Vec<Scalar>) -> DenseRowMatrix<Scalar> {
    let x = v[0];
    let y = v[1];
    DenseRowMatrix::new(2, vec![-1.0, (y + 0.5).cos(), -(x - 2.0).sin(), 1.0])
}

fn test_newton_method_n() {
    println!("\x1b[32m=== Vector newton method ===\x1b[0m");
    let initial_x = vec![0.0, 0.0];
    let max_iterations = 100;
    let answer = newton_method_n(task_fn, task_fnprime, initial_x, 1e-4, max_iterations);
    match answer {
        Some(x) => {
            println!("Result = {x:?}");
            println!("Residual = {r:?}", r = task_fn(&x));
        }
        None => println!("Newton method did not terminate in {max_iterations}"),
    }
}

fn task_fi(lambda: Scalar, v: &Vec<Scalar>) -> Vec<Scalar> {
    let x = v[0];
    let y = v[1];
    vec![
        -x + lambda * (y + 0.5).sin() - 1.0,
        lambda * (x - 2.0).cos() + y,
    ]
}

fn task_fiprime(lambda: Scalar, v: &Vec<Scalar>) -> DenseRowMatrix<Scalar> {
    let x = v[0];
    let y = v[1];
    DenseRowMatrix::new(
        2,
        vec![
            -1.0,
            lambda * (y + 0.5).cos(),
            -lambda * (x - 2.0).sin(),
            1.0,
        ],
    )
}

fn test_newton_method_interpolated() {
    println!("\x1b[32m=== Vector interpolated newton method ===\x1b[0m");
    let max_iterations = 100;
    let answer = newton_method_n_interpolated(task_fi, task_fiprime, 2, 10, 1e-4, max_iterations);
    match answer {
        Some(x) => {
            println!("Result = {x:?}");
            println!("Residual = {r:?}", r = task_fn(&x));
        }
        None => println!("Newton method did not terminate in {max_iterations}"),
    }
}

fn main() {
    test_newton_method();
    test_newton_method_n();
    test_newton_method_interpolated();
}
