use basic::{IntInterval, Interval, Scalar};
use localize_root::localize_root_in_interval;
use newton_method::newton_method;

pub mod basic;
pub mod localize_root;
pub mod newton_method;

fn task_f(x: Scalar) -> Scalar {
    (Scalar::from(0.5) * x + Scalar::from(0.2)).tan() - x * x
}

fn task_fprime(x: Scalar) -> Scalar {
    Scalar::from(1.0)
        / (Scalar::from(2.0) * (x * Scalar::from(0.5) + Scalar::from(0.2)).cos().powi(2))
        - Scalar::from(2.0) * x
}

fn main() {
    // let interval = localize_root_in_interval(task_f, IntInterval::new(0, 16));
    newton_method(task_f, task_fprime, Interval::new(0.0, 9.0), 1e-10, 100);
    // println!(
    //     "interval={:?}",
    //     localize_root(|x| 10.0 * (x - 0.6) * (x - 0.6) - 1.0)
    // );
}
