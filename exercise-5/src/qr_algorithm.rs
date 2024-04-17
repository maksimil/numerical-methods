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

pub fn hessenberg_form(matrix: &mut Vec<Scalar>, dimension: Index) {
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

        if length == SCALAR_ZERO {
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
