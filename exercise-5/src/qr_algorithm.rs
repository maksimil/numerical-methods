use crate::scalar::*;

pub fn compute_reflection(
    vector: &[Scalar],
    direction: &mut [Scalar],
    scale: &mut Scalar,
) -> Scalar {
    let dimension = vector.len();

    let length = {
        let mut r = SCALAR_ZERO;

        for k in 0..dimension {
            r += vector[k] * vector[k];
        }

        r.sqrt()
    };

    direction[0] = vector[0] - length;

    for k in 1..dimension {
        direction[k] = vector[k];
    }

    *scale = SCALAR_ZERO;

    for k in 0..dimension {
        *scale += direction[k] * direction[k];
    }

    length
}

pub fn hessenberg_form(matrix: &mut [Scalar], dimension: Index) {
    if dimension < 3 {
        return;
    }

    let mut direction = vec![SCALAR_ZERO; dimension];
    let mut scale = SCALAR_ZERO;

    let mut working_column = vec![SCALAR_ZERO; dimension];

    for col in 0..dimension - 2 {
        // computing the reflection
        for k in col + 1..dimension {
            working_column[k] = matrix[dimension * k + col];
        }

        let length = compute_reflection(
            &working_column[col + 1..dimension],
            &mut direction[col + 1..dimension],
            &mut scale,
        );

        if scale == SCALAR_ZERO {
            continue;
        }

        // applying reflection to the matrix
        // A <- HA
        matrix[dimension * (col + 1) + col] = length;

        for r in col + 2..dimension {
            matrix[dimension * r + col] = SCALAR_ZERO;
        }

        for c in col + 1..dimension {
            let coeficient = {
                let mut res = SCALAR_ZERO;

                for k in col + 1..dimension {
                    res += direction[k] * matrix[dimension * k + c];
                }

                -(2 as Scalar) * res / scale
            };

            for r in col + 1..dimension {
                matrix[dimension * r + c] += coeficient * direction[r];
            }
        }

        // A <- AH
        for r in 0..dimension {
            let coeficient = {
                let mut res = SCALAR_ZERO;

                for k in col + 1..dimension {
                    res += direction[k] * matrix[dimension * r + k];
                }

                -(2 as Scalar) * res / scale
            };

            for c in col + 1..dimension {
                matrix[dimension * r + c] += coeficient * direction[c];
            }
        }
    }
}

// A=QR -> RQ
// A is in hessenberg form
pub fn qr_iteration(matrix: &mut [Scalar], dimension: Index) {
    let mut reflections = vec![SCALAR_ZERO; 2 * (dimension - 1)];
    let mut scales = vec![SCALAR_ZERO; dimension - 1];
    let mut work_vector = vec![SCALAR_ZERO; 2];

    // computing QR, A <- R
    for col in 0..dimension - 1 {
        // reflecting
        work_vector[0] = matrix[dimension * col + col];
        work_vector[1] = matrix[dimension * (col + 1) + col];

        let length = compute_reflection(
            &work_vector,
            &mut reflections[col * 2..col * 2 + 2],
            &mut scales[col],
        );

        let scale = scales[col];
        let reflection = &reflections[col * 2..col * 2 + 2];

        if scale == SCALAR_ZERO {
            continue;
        }

        // applying the reflection
        matrix[dimension * col + col] = length;
        matrix[dimension * (col + 1) + col] = SCALAR_ZERO;

        for c in col + 1..dimension {
            let coeficient = -(2 as Scalar)
                * (reflection[0] * matrix[dimension * col + c]
                    + reflection[1] * matrix[dimension * (col + 1) + c])
                / scale;

            matrix[dimension * col + c] += coeficient * reflection[0];
            matrix[dimension * (col + 1) + c] += coeficient * reflection[1];
        }
    }

    // println!("R = {matrix:#?}");

    // computing RQ, A <- AQ
    for k in 0..dimension - 1 {
        let scale = scales[k];
        let reflection = &reflections[k * 2..k * 2 + 2];

        if scale == SCALAR_ZERO {
            continue;
        }

        for r in 0..dimension {
            let coeficient = -(2 as Scalar)
                * (reflection[0] * matrix[dimension * r + k]
                    + reflection[1] * matrix[dimension * r + k + 1])
                / scale;

            matrix[dimension * r + k] += coeficient * reflection[0];
            matrix[dimension * r + k + 1] += coeficient * reflection[1];
        }
    }
}

const SUBDIAGONAL_SMALL: Scalar = 1e-8;

pub fn qr_algorithm(matrix: &[Scalar], dimension: Index) -> Option<Vec<Scalar>> {
    let mut work_matrix = matrix.to_owned();
    let mut work_dimension = dimension;
    let mut values = vec![];
    values.reserve(dimension);

    hessenberg_form(&mut work_matrix, dimension);

    let mut did_converge = false;

    for _ in 0..MAX_ITERATIONS {
        // iterations
        let shift = work_matrix[work_dimension * work_dimension - 1];

        for k in 0..work_dimension {
            work_matrix[work_dimension * k + k] -= shift;
        }

        qr_iteration(&mut work_matrix, work_dimension);

        for k in 0..work_dimension {
            work_matrix[work_dimension * k + k] += shift;
        }

        // dimension reduction
        if work_matrix[work_dimension * work_dimension - 2].abs() < SUBDIAGONAL_SMALL {
            values.push(work_matrix[work_dimension * work_dimension - 1]);

            if work_dimension == 2 {
                values.push(work_matrix[0]);
                did_converge = true;
                break;
            }

            for r in 1..work_dimension - 1 {
                for c in 0..work_dimension - 1 {
                    work_matrix[(work_dimension - 1) * r + c] = work_matrix[work_dimension * r + c];
                }
            }

            work_dimension -= 1;

            work_matrix.resize(work_dimension * work_dimension, SCALAR_ZERO);
        }
    }

    values.sort_by(Scalar::total_cmp);

    if did_converge {
        Some(values)
    } else {
        None
    }
}
