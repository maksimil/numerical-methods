use crate::{
    basic::{Interval, Scalar},
    localize_root::check_root,
};

fn newton_iteration(
    f: &impl Fn(Scalar) -> Scalar,
    fprime: &impl Fn(Scalar) -> Scalar,
    x: Scalar,
) -> Scalar {
    let fx = f(x);
    let fprimex = fprime(x);
    let x_new = x - fx / fprimex;
    println!("f(x)={fx}, fprime(x)={fprimex}, x={x}, xnew={x_new}");
    x - f(x) / fprime(x)
}

#[derive(PartialEq, Eq)]
enum SplitAction {
    None,
    ShrinkToStart,
    ShrinkToEnd,
}

fn interval_split(
    f: &impl Fn(Scalar) -> Scalar,
    interval: Interval,
    x: Scalar,
) -> (Interval, SplitAction) {
    let left_interval = Interval::new(interval.start, x);
    let right_interval = Interval::new(x, interval.end);
    if check_root(f, &left_interval) {
        (left_interval, SplitAction::ShrinkToStart)
    } else {
        (right_interval, SplitAction::ShrinkToEnd)
    }
}

pub fn newton_method(
    f: impl Fn(Scalar) -> Scalar,
    fprime: impl Fn(Scalar) -> Scalar,
    mut interval: Interval,
    accuracy: Scalar,
    max_iterations: usize,
) -> Option<Scalar> {
    let mut previous_iteration = SplitAction::None;

    for k in 0..max_iterations {
        println!("k={k}");
        // from interval end
        if previous_iteration != SplitAction::ShrinkToEnd {
            let x = interval.end;
            let x_new = newton_iteration(&f, &fprime, x);

            if interval.contains(x_new) {
                if (x_new - x).abs() < accuracy {
                    return Some(x_new);
                } else {
                    (interval, previous_iteration) = interval_split(&f, interval, x_new);
                    println!("New interval={:?}", interval);
                    continue;
                }
            }
        }

        // from interval start
        if previous_iteration != SplitAction::ShrinkToStart {
            let x = interval.start;
            let x_new = newton_iteration(&f, &fprime, x);

            if interval.contains(x_new) {
                if (x_new - x).abs() < accuracy {
                    return Some(x_new);
                } else {
                    (interval, previous_iteration) = interval_split(&f, interval, x_new);
                    println!("New interval={:?}", interval);
                    continue;
                }
            }
        }

        // split interval
        {
            println!("Overshot");
            let x_new = interval.middle();
            (interval, previous_iteration) = interval_split(&f, interval, x_new);
            println!("New interval={:?}", interval);
            if interval.end - interval.start < accuracy {
                return Some(interval.middle());
            }
            continue;
        }
    }

    None
}
