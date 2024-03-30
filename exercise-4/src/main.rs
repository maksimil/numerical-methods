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
    let graph_m = 10_000;
    let ns = vec![3, 5, 10, 15, 20, 25, 40, 50, 60, 70, 80, 90, 100];
    let graph_ni = vec![0, 1, 2, 4, 7, 12];

    let mut stats_file = File::create("output/stats.csv").unwrap();
    let mut graph_file = File::create("output/graph.csv").unwrap();

    let mut lagranges = vec![];
    let mut lagranges_opt = vec![];
    let mut newtons = vec![];
    let mut newtons_opt = vec![];
    let mut spline1s = vec![];
    let mut spline2s = vec![];
    let mut spline3s = vec![];

    stats_file
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

        let spline1 = Spline::new(1, &points, example_function);
        let spline2 = Spline::new(2, &points, example_function);
        let spline3 = Spline::new(3, &points, example_function);

        let lagrange_err = calc_error(example_function, |x| lagrange_poly.at(x), a, b, m);
        let lagrange_opt_err = calc_error(example_function, |x| lagrange_opt_poly.at(x), a, b, m);

        let newton_err = calc_error(example_function, |x| newton_poly.at(x), a, b, m);
        let newton_opt_err = calc_error(example_function, |x| newton_opt_poly.at(x), a, b, m);

        let spline1_err = calc_error(example_function, |x| spline1.at(x), a, b, m);
        let spline2_err = calc_error(example_function, |x| spline2.at(x), a, b, m);
        let spline3_err = calc_error(example_function, |x| spline3.at(x), a, b, m);

        stats_file
            .write(
                format!(
                    "{n};{m};{lagrange_err:.e};{lagrange_opt_err:.e};{newton_err:.e};\
                        {newton_opt_err:.e};{spline1_err:.e};{spline2_err:.e};{spline3_err:.e}\n"
                )
                .as_bytes(),
            )
            .unwrap();

        lagranges.push(lagrange_poly);
        lagranges_opt.push(lagrange_opt_poly);

        newtons.push(newton_poly);
        newtons_opt.push(newton_opt_poly);

        spline1s.push(spline1);
        spline2s.push(spline2);
        spline3s.push(spline3);
    }

    graph_file.write("x;f".as_bytes()).unwrap();
    for ii in 0..graph_ni.len() {
        let i = graph_ni[ii];
        let n = ns[i];
        graph_file
            .write(format!(";L_{n};RL_{n};Lopt_{n};RLopt_{n};N_{n};RN_{n};Nopt_{n};RNopt_{n};S1_{n};RS1_{n};S2_{n};RS2_{n};S3_{n};RS3_{n}").as_bytes())
            .unwrap();
    }
    graph_file.write("\n".as_bytes()).unwrap();

    for k in 0..=graph_m {
        let x = a + (b - a) * (k as Scalar) / (graph_m as Scalar);
        let f = example_function(x);
        graph_file.write(format!("{x};{f}").as_bytes()).unwrap();

        for ii in 0..graph_ni.len() {
            let i = graph_ni[ii];
            let l = lagranges[i].at(x);
            let lopt = lagranges_opt[i].at(x);
            let nt = newtons[i].at(x);
            let nt_opt = newtons_opt[i].at(x);
            let s1 = spline1s[i].at(x);
            let s2 = spline2s[i].at(x);
            let s3 = spline3s[i].at(x);
            graph_file
                .write(format!(";{l};{dl};{lopt};{dlopt};{nt};{dnt};{nt_opt};{dnt_opt};{s1};{ds1};{s2};{ds2};{s3};{ds3}",
                               dl=f-l, dlopt = f-lopt, dnt=f-nt, dnt_opt = f-nt_opt, ds1=f-s1,ds2=f-s2,ds3=f-s3).as_bytes())
                .unwrap();
        }

        graph_file.write("\n".as_bytes()).unwrap();
    }
}
