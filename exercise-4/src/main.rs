use std::{fs::File, io::Write};

use spline::Spline;

use crate::{
    points::{points_cheb, points_eq},
    polynomial::{create_lagrange, create_newton},
    scalar::*,
};

mod lu;
mod points;
mod polynomial;
mod scalar;
mod spline;

fn example_function(x: Scalar) -> Scalar {
    (0.5 * x + 0.2).tan() - x * x
}

fn calc_error(
    f: impl Fn(Scalar) -> Scalar,
    g: impl Fn(Scalar) -> Scalar,
    a: Scalar,
    b: Scalar,
    m: Index,
) -> Scalar {
    let mut r = SCALAR_ZERO;
    for k in 0..=m {
        let x = a + (b - a) * (k as Scalar) / (m as Scalar);
        let d = (f(x) - g(x)).abs();
        if d > r {
            r = d;
        }
    }
    r
}

fn main() {
    let a = -1.;
    let b = 1.;

    let m = 10_000;
    let ns = vec![3, 5, 10, 15, 20, 25, 40, 50, 60, 70, 80, 90, 100];

    let initial = vec![2.548, -2.169];

    let mut poly_stats_file = File::create("output/stats.csv").unwrap();

    poly_stats_file
        .write("n;m;RL;RLopt;RN;RNopt;RS1;RS2;RS3\n".as_bytes())
        .unwrap();

    for i in 0..ns.len() {
        let n = ns[i];

        let points = points_eq(a, b, n);
        let points_opt = points_cheb(a, b, n);

        let lagrange_poly = create_lagrange(&points, example_function);
        let lagrange_opt_poly = create_lagrange(&points_opt, example_function);

        let newton_poly = create_newton(&points, example_function);
        let newton_opt_poly = create_newton(&points_opt, example_function);

        let spline1 = Spline::new(1, &points, example_function, &initial);
        let spline2 = Spline::new(2, &points, example_function, &initial);
        let spline3 = Spline::new(3, &points, example_function, &initial);

        let lagrange_err = calc_error(example_function, |x| lagrange_poly.at(x), a, b, m);
        let lagrange_opt_err = calc_error(example_function, |x| lagrange_opt_poly.at(x), a, b, m);

        let newton_err = calc_error(example_function, |x| newton_poly.at(x), a, b, m);
        let newton_opt_err = calc_error(example_function, |x| newton_opt_poly.at(x), a, b, m);

        let spline1_err = calc_error(example_function, |x| spline1.at(x), a, b, m);
        let spline2_err = calc_error(example_function, |x| spline2.at(x), a, b, m);
        let spline3_err = calc_error(example_function, |x| spline3.at(x), a, b, m);

        poly_stats_file
            .write(
                format!(
                    "{n};{m};{lagrange_err:.e};{lagrange_opt_err:.e};{newton_err:.e};\
                        {newton_opt_err:.e};{spline1_err:.e};{spline2_err:.e};{spline3_err:.e}\n"
                )
                .as_bytes(),
            )
            .unwrap();
    }
}
