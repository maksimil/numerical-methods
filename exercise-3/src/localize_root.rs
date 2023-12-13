use crate::basic::{IntInterval, Interval, Scalar};

pub fn check_root(f: impl Fn(Scalar) -> Scalar, interval: &Interval) -> bool {
    let start_value = f(interval.start);
    let end_value = f(interval.end);
    return start_value * end_value <= Scalar::from(0);
}

#[allow(dead_code)]
fn log_parameters<IF, PF>(k_max: u32, interval_fn: &IF, parts_fn: &PF)
where
    IF: Fn(u32) -> IntInterval,
    PF: Fn(u32) -> u32,
{
    println!(
        "Max checked intervals = {:.e}",
        (1..=k_max)
            .map(|k| {
                let interval = interval_fn(k);
                u64::from(interval.diameter()) * u64::from(parts_fn(k))
            })
            .sum::<u64>(),
    );

    for k in 1..=k_max {
        let interval = interval_fn(k);
        let interval_parts = parts_fn(k);
        let delta = Scalar::from(1) / Scalar::from(interval_parts);
        let count = interval.diameter() * interval_parts;
        println!("k={k}, interval=[{interval_start:.e}, {interval_end:.e}], parts={interval_parts:.e}, delta={delta:.e} count={count:.e}", 
                 interval_start=interval.start, interval_end = interval.end);
    }
}

pub fn general_localize_root(
    f: impl Fn(Scalar) -> Scalar + Clone,
    k_max: u32,
    interval_fn: impl Fn(u32) -> IntInterval,
    parts_fn: impl Fn(u32) -> u32,
) -> Option<Interval> {
    let mut k = 1;

    // log_parameters(k_max, &interval_fn, &parts_fn);

    while k <= k_max {
        let general_interval = interval_fn(k);
        let interval_parts = parts_fn(k);

        for i in 0..general_interval.diameter() {
            for j in 0..interval_parts {
                let current_interval = Interval::new(
                    Scalar::from(general_interval.start)
                        + Scalar::from(i)
                        + Scalar::from(j) / Scalar::from(interval_parts),
                    Scalar::from(general_interval.start)
                        + Scalar::from(i)
                        + Scalar::from(j + 1) / Scalar::from(interval_parts),
                );
                // println!("Checking {current_interval:?}");

                if check_root(f.clone(), &current_interval) {
                    return Some(current_interval);
                }
            }
        }

        k += 1;
    }

    None
}

// Localizes roots in radius=100_000 and dx=3e-2 in about 200ms
// for a const function
pub fn localize_root(f: impl Fn(Scalar) -> Scalar + Clone) -> Option<Interval> {
    general_localize_root(
        f,
        3,
        |k| {
            let radius = 10i32.pow(k - 1);
            IntInterval::new(-radius, radius)
        },
        |k| 2u32.pow(k + 2),
    )
}

pub fn localize_root_in_interval(
    f: impl Fn(Scalar) -> Scalar + Clone,
    interval: IntInterval,
) -> Option<Interval> {
    general_localize_root(f, 3, |_| interval.clone(), |k| 100u32.pow(k))
}
