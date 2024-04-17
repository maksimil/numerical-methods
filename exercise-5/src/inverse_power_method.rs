use crate::{
    lu::LUDecomposition,
    power_method::{PowerMethod, PM_CONVERGENCE, PM_MIN_COORD},
    scalar::*,
};

pub fn inverse_power_method(
    matrix: &[Scalar],
    dimension: Index,
    initial: Scalar,
) -> Option<PowerMethod> {
    let mut y_vector = vec![SCALAR_ZERO; dimension];
    let mut z_vector = vec![SCALAR_ZERO; dimension];
    let mut system_matrix = vec![SCALAR_ZERO; dimension * dimension];
    let mut shift = initial;
    let mut stop_shifting = false;
    let mut lu = LUDecomposition::new_zero(0);

    y_vector[0] = SCALAR_ONE;

    for i in 0..dimension {
        for j in 0..dimension {
            system_matrix[dimension * i + j] = matrix[dimension * i + j];
        }
    }

    let mut did_converge = false;

    for k in 0..MAX_ITERATIONS {
        println!("k={k}, shift={shift}");
        // setting z=y (y is normed)
        for k in 0..dimension {
            z_vector[k] = y_vector[k];
        }

        // computing, (A-sE)y=z
        if !stop_shifting {
            for k in 0..dimension {
                system_matrix[dimension * k + k] = matrix[dimension * k + k] - shift;
            }

            let lu_opt = LUDecomposition::compute(&system_matrix, dimension);

            match lu_opt {
                Some(mut lu_computed) => {
                    std::mem::swap(&mut lu, &mut lu_computed);
                }
                None => {
                    println!("Stopped shifting");
                    stop_shifting = true;
                }
            };
        }

        lu.solve(&mut y_vector);

        // computing mu
        let mu = if !stop_shifting {
            let mut count = 0;
            let mut mu = SCALAR_ZERO;

            for k in 0..dimension {
                if y_vector[k].abs() >= PM_MIN_COORD {
                    count += 1;
                    mu += z_vector[k] / y_vector[k];
                }
            }

            mu / (count as Scalar)
        } else {
            SCALAR_ZERO
        };

        // norming, y = ort(y)
        let norm = {
            let mut r = SCALAR_ZERO;

            for k in 0..dimension {
                r += y_vector[k].abs();
            }

            r
        };

        for k in 0..dimension {
            y_vector[k] = y_vector[k] / norm;
        }

        // checking convergence of shift and z (z next is y)
        let z_norm = {
            let mut r = SCALAR_ZERO;

            for k in 0..dimension {
                r += z_vector[k].abs()
            }

            r
        };

        let y_norm = {
            let mut r = SCALAR_ZERO;

            for k in 0..dimension {
                r += y_vector[k].abs();
            }

            r
        };

        let diff_norm = {
            let mut r = SCALAR_ZERO;

            let y_sign = if y_vector[0] * z_vector[0] > SCALAR_ZERO {
                SCALAR_ONE
            } else {
                -SCALAR_ONE
            };

            for k in 0..dimension {
                r += (y_sign * y_vector[k] - z_vector[k]).abs();
            }

            r
        };

        if diff_norm < PM_CONVERGENCE * z_norm.max(y_norm)
            && (mu < PM_CONVERGENCE * shift.abs().max((shift + mu).abs()) || stop_shifting)
        {
            did_converge = true;
            break;
        } else {
            shift += mu;
        }
    }

    if did_converge {
        Some(PowerMethod {
            value: shift,
            vector: y_vector,
        })
    } else {
        None
    }
}
