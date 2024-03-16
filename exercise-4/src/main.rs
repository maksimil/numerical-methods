use polynomial::Scalar;

use crate::polynomial::{create_lagrange, create_newton};

mod polynomial;

fn example_function(x: Scalar) -> Scalar {
    return x * x * x + 1.0;
}

fn main() {
    let points = vec![-1.0, 0.0, 1.0, 10.0, 3.0];

    let lagrange = create_lagrange(&points, example_function);
    println!("Lagrange: {lagrange:?}");

    let newton = create_newton(&points, example_function);
    println!("Newton: {newton:?}");
}
