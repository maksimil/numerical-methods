use crate::{lu::LUDecomposition, scalar::*};
use rand;

#[derive(Debug, Clone)]
pub struct GenerateMatrix {
    pub matrix: Vec<Scalar>,
    pub eigenvalues: Vec<Scalar>,
}

const VALUES_SEPARATION: Scalar = 0.5;
const VALUES_MAGNITUDE_MUL: Scalar = 10.;
const TRANSFROMATION_VALUES_MAGNITUDE: Scalar = 5.;

fn rand_scalar(magnitude: Scalar) -> Scalar {
    (rand::random::<Scalar>() - 0.5) * 2. * magnitude
}

pub fn generate_matrix(dimension: Index) -> GenerateMatrix {
    let values_between = VALUES_SEPARATION;
    let values_magnitude = VALUES_MAGNITUDE_MUL * values_between * (dimension as Scalar);
    let transformation_matrix_values_magnitude = TRANSFROMATION_VALUES_MAGNITUDE;

    let eigenvalues = {
        let mut values = vec![SCALAR_ZERO; dimension];

        for k in 0..dimension {
            let mut singular_value = SCALAR_ZERO;
            let mut valid = false;

            while !valid {
                singular_value = rand_scalar(values_magnitude);

                valid = true;
                for l in 0..k {
                    if (singular_value - values[l]).abs() < values_between {
                        valid = false;
                        break;
                    }
                }
            }

            values[k] = singular_value;
        }

        values
    };

    let (transformation, transformation_lu) = {
        let mut matrix = vec![SCALAR_ZERO; dimension * dimension];
        let mut lu = None;

        while lu.is_none() {
            for k in 0..dimension * dimension {
                matrix[k] = rand_scalar(transformation_matrix_values_magnitude);
            }

            println!("Checking {:?}", matrix);

            lu = LUDecomposition::compute(&matrix, dimension);
        }

        (matrix, lu.unwrap())
    };

    let matrix = {
        let mut matrix = vec![SCALAR_ZERO; dimension * dimension];
        let mut work_array = vec![SCALAR_ZERO; dimension];

        for col in 0..dimension {
            for row in 0..dimension {
                work_array[row] = eigenvalues[row] * transformation[row * dimension + col];
            }

            transformation_lu.solve(&mut work_array);

            for row in 0..dimension {
                matrix[row * dimension + col] = work_array[row];
            }
        }

        matrix
    };

    GenerateMatrix {
        matrix,
        eigenvalues,
    }
}
