use crate::scalar::*;
use spline::Spline;

mod lu;
mod polynomial;
mod scalar;
mod spline;

fn example_function(x: Scalar) -> Scalar {
    return x * x + 1.0;
}

fn main() {
    let s = Spline::new(1, vec![-2., -1., 0., 1., 2.], example_function);
    println!("{s:?}");
    // let points = vec![-1.0, 0.0, 1.0, 10.0, 3.0];
    //
    // let lagrange = create_lagrange(&points, example_function);
    // println!("Lagrange: {lagrange:?}");
    //
    // let newton = create_newton(&points, example_function);
    // println!("Newton: {newton:?}");
}
