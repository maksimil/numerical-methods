use exercise_2::{
    lu_decomposition::LUDecomposition,
    matrix::{
        column::{ColumnFuncInitializer, ColumnMut, ColumnRef},
        dense::DenseRowMatrix,
        norms::NormedColumn,
    },
};

use crate::basic::Scalar;

// returns -\Delta x
fn newton_diff(
    f: &impl Fn(&Vec<Scalar>) -> Vec<Scalar>,
    fprime: &impl Fn(&Vec<Scalar>) -> DenseRowMatrix<Scalar>,
    x: &Vec<Scalar>,
) -> Vec<Scalar> {
    let jacobian = fprime(x);
    let decomposition = LUDecomposition::calculate(jacobian);
    let mut v = f(x);
    decomposition.solve(&mut v);
    v
}

pub fn newton_method_n(
    f: impl Fn(&Vec<Scalar>) -> Vec<Scalar>,
    fprime: impl Fn(&Vec<Scalar>) -> DenseRowMatrix<Scalar>,
    mut x: Vec<Scalar>,
    accuracy: Scalar,
    max_steps: usize,
) -> Option<Vec<Scalar>> {
    for k in 0..max_steps {
        println!("k={k}");
        let x_diff = newton_diff(&f, &fprime, &x);
        let x_new = Vec::<Scalar>::new_func(x.dimension(), |i| x.at(i) - x_diff.at(i));

        println!("New x = {x_new:?}, residue {r:?}", r = f(&x_new));

        if x_diff.norm_inf() < accuracy {
            return Some(x_new);
        } else {
            x = x_new;
        }
    }

    None
}

pub fn newton_method_n_staged(
    f: impl Fn(u32, &Vec<Scalar>) -> Vec<Scalar>,
    fprime: impl Fn(u32, &Vec<Scalar>) -> DenseRowMatrix<Scalar>,
    dimension: usize,
    stages: u32,
    accuracy: Scalar,
    max_steps: usize,
) -> Option<Vec<Scalar>> {
    println!("stage = 0");
    let mut x = {
        let mut initial_dx = newton_diff(
            &|v| f(0, v),
            &|v| fprime(0, v),
            &Vec::new_fill(dimension, 0.0),
        );
        for k in 0..dimension {
            *initial_dx.at_mut(k) = -initial_dx.at(k);
        }
        initial_dx
    };
    println!("Initial estimate = {x:?}");

    for stage in 1..=stages {
        println!("stage = {stage}");

        let x_new_result = newton_method_n(
            |v| f(stage, v),
            |v| fprime(stage, v),
            x.clone(),
            accuracy,
            max_steps,
        );

        println!("New value = {x_new_result:?}");

        match x_new_result {
            Some(x_new) => {
                x = x_new;
            }
            None => return None,
        }
    }

    Some(x)
}

pub fn newton_method_n_interpolated(
    f: impl Fn(Scalar, &Vec<Scalar>) -> Vec<Scalar>,
    fprime: impl Fn(Scalar, &Vec<Scalar>) -> DenseRowMatrix<Scalar>,
    dimension: usize,
    stages: u32,
    accuracy: Scalar,
    max_steps: usize,
) -> Option<Vec<Scalar>> {
    newton_method_n_staged(
        |stage, v| f(Scalar::from(stage) / Scalar::from(stages), v),
        |stage, v| fprime(Scalar::from(stage) / Scalar::from(stages), v),
        dimension,
        stages,
        accuracy,
        max_steps,
    )
}
