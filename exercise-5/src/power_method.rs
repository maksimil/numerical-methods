use crate::{lu::dot_product, scalar::*};

#[derive(Debug, Clone)]
pub struct PowerMethod {
    pub value: Scalar,
    pub vector: Vec<Scalar>,
}

const PM_MIN_COORD: Scalar = 1e-8;
const PM_CONVERGENCE: Scalar = 1e-6;

pub fn power_method(matrix: &[Scalar], dimension: Index) -> Option<PowerMethod> {
    let mut y_vector = vec![SCALAR_ZERO; dimension];
    let mut z_vector = vec![SCALAR_ZERO; dimension];
    y_vector[0] = SCALAR_ONE;

    let mut previous_lambda = vec![SCALAR_ZERO; dimension];
    let mut previous_lambda_valid = vec![false; dimension];
    let mut current_lambda = vec![SCALAR_ZERO; dimension];
    let mut current_lambda_valid = vec![false; dimension];

    let mut did_converge = false;

    for _ in 0..MAX_ITERATIONS {
        // norming the vector
        let norm = {
            let mut r = SCALAR_ZERO;

            for i in 0..dimension {
                r += y_vector[i].abs();
            }

            r
        };

        for i in 0..dimension {
            z_vector[i] = y_vector[i] / norm;
        }

        // computing y_vector, y=Az
        for i in 0..dimension {
            y_vector[i] = dot_product(&matrix[i * dimension..i * dimension + dimension], &z_vector);
        }

        // computing lambdas
        for i in 0..dimension {
            previous_lambda[i] = current_lambda[i];
            previous_lambda_valid[i] = current_lambda_valid[i];

            if z_vector[i].abs() < PM_MIN_COORD {
                current_lambda_valid[i] = false;
                current_lambda[i] = SCALAR_ZERO;
            } else {
                current_lambda_valid[i] = true;
                current_lambda[i] = y_vector[i] / z_vector[i];
            }
        }

        // checking convergence
        let prev_lambda_norm = {
            let mut r = SCALAR_ZERO;

            for i in 0..dimension {
                if previous_lambda_valid[i] {
                    r += previous_lambda[i].abs();
                }
            }

            r
        };

        let cur_lambda_norm = {
            let mut r = SCALAR_ZERO;

            for i in 0..dimension {
                if current_lambda_valid[i] {
                    r += current_lambda[i].abs();
                }
            }

            r
        };

        let diff_lambda_norm = {
            let mut r = SCALAR_ZERO;

            for i in 0..dimension {
                if previous_lambda_valid[i] || current_lambda_valid[i] {
                    r += (previous_lambda[i] - current_lambda[i]).abs();
                }
            }

            r
        };

        if diff_lambda_norm < PM_CONVERGENCE * prev_lambda_norm.max(cur_lambda_norm) {
            did_converge = true;
            break;
        }
    }

    if did_converge {
        // compute value
        let value = {
            let mut r = SCALAR_ZERO;
            let mut k = 0;

            for i in 0..dimension {
                if current_lambda_valid[i] {
                    k += 1;
                    r += current_lambda[i];
                }
            }

            r / (k as Scalar)
        };

        // norming the vector
        let norm = {
            let mut r = SCALAR_ZERO;

            for i in 0..dimension {
                r += y_vector[i].abs();
            }

            r
        };

        for i in 0..dimension {
            z_vector[i] = y_vector[i] / norm;
        }

        Some(PowerMethod {
            value,
            vector: z_vector,
        })
    } else {
        None
    }
}
