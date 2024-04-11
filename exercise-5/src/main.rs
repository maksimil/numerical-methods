use crate::{lu::LUDecomposition, scalar::*};

mod lu;
mod scalar;

fn main() {
    let dimension = 2;
    let values_between = 0.5;
    let values_magnitude = ((10 * dimension) as Scalar) * values_between;
    let transformation_matrix_values_magnitude = 5.;

    let diagonal = {
        let mut diagonal = vec![SCALAR_ZERO; dimension];

        for k in 0..dimension {
            let mut singular_value = SCALAR_ZERO;
            let mut valid = false;

            while !valid {
                singular_value = (rand::random::<Scalar>() - 0.5) * 2. * values_magnitude;

                valid = true;
                for l in 0..k {
                    if (singular_value - diagonal[l]).abs() < values_between {
                        valid = false;
                        break;
                    }
                }
            }

            diagonal[k] = singular_value;
        }

        diagonal
    };

    let (transformation, transformation_lu) = {
        let mut matrix = vec![SCALAR_ZERO; dimension * dimension];
        let mut lu = None;

        while lu.is_none() {
            for k in 0..dimension * dimension {
                matrix[k] =
                    (rand::random::<Scalar>() - 0.5) * 2. * transformation_matrix_values_magnitude;
            }

            println!("Checking {:?}", matrix);

            lu = LUDecomposition::compute(&matrix, dimension);
        }

        (matrix, lu.unwrap())
    };

    println!("{:?}", diagonal);
    println!("{:?}", transformation);

    let matrix = {
        let mut matrix = vec![SCALAR_ZERO; dimension * dimension];
        let mut work_array = vec![SCALAR_ZERO; dimension];

        for col in 0..dimension {
            for row in 0..dimension {
                work_array[row] = diagonal[row] * transformation[row * dimension + col];
            }

            transformation_lu.solve(&mut work_array);

            for row in 0..dimension {
                matrix[row * dimension + col] = work_array[row];
            }
        }

        matrix
    };

    println!("{:?}", matrix);
}
