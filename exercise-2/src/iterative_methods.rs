use crate::{
    basic::{Numerical, OtherNumericalOps},
    matrix::{
        column::{apply, ColumnFunc, ColumnFuncInitializer, ColumnMut, ColumnRef},
        dense::DenseRowMatrix,
        norms::NormedMatrix,
        traits::{MatrixFunc, MatrixFuncInitializer, MatrixRef},
        transpose::MatrixSymmetricSquare,
    },
};

pub fn iterative_solve<Matrix, ColumnIn, ColumnOut>(
    matrix: &Matrix,
    vector: &ColumnIn,
    accuracy: Matrix::Scalar,
) -> ColumnOut
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnMut<Scalar = Matrix::Scalar> + ColumnFuncInitializer<Scalar = Matrix::Scalar>,
{
    internal_iterative_solve(matrix, vector, accuracy, false)
}

fn internal_iterative_solve<Matrix, ColumnIn, ColumnOut>(
    a_matrix: &Matrix,
    vector: &ColumnIn,
    accuracy: Matrix::Scalar,
    can_abuse_norm: bool,
) -> ColumnOut
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnMut<Scalar = Matrix::Scalar> + ColumnFuncInitializer<Scalar = Matrix::Scalar>,
{
    let dimension = a_matrix.dimension();
    let e_matrix = MatrixFunc::new(dimension, |row, column| {
        if row == column {
            Matrix::Scalar::one()
        } else {
            Matrix::Scalar::zero()
        }
    });

    // try p=1
    {
        let a_norm_one = a_matrix.norm_one();
        let b_matrix = MatrixFunc::new(dimension, |row, column| {
            e_matrix.at(row, column) - a_matrix.at(row, column) / a_norm_one
        });
        let q_one = b_matrix.norm_one();

        if q_one < Matrix::Scalar::one() {
            let b_matrix: DenseRowMatrix<Matrix::Scalar> = DenseRowMatrix::from_matrix(&b_matrix);
            let c_vector: Vec<Matrix::Scalar> =
                Vec::from_column(&ColumnFunc::new(dimension, |i| vector.at(i) / a_norm_one));
            return run_method(b_matrix, NormEnum::One { value: q_one }, c_vector, accuracy);
        }
    }

    // try p=inf
    {
        let a_norm_inf = a_matrix.norm_inf();
        let b_matrix = MatrixFunc::new(dimension, |row, column| {
            e_matrix.at(row, column) - a_matrix.at(row, column) / a_norm_inf
        });
        let q_inf = b_matrix.norm_inf();

        if q_inf < Matrix::Scalar::one() {
            let b_matrix: DenseRowMatrix<Matrix::Scalar> = DenseRowMatrix::from_matrix(&b_matrix);
            let c_vector: Vec<Matrix::Scalar> =
                Vec::from_column(&ColumnFunc::new(dimension, |i| vector.at(i) / a_norm_inf));
            return run_method(
                b_matrix,
                NormEnum::Infty { value: q_inf },
                c_vector,
                accuracy,
            );
        }
    }

    // abuse norm if can
    if can_abuse_norm {
        let a_norm_one = a_matrix.norm_one();
        let b_matrix: DenseRowMatrix<Matrix::Scalar> =
            DenseRowMatrix::from_matrix(&MatrixFunc::new(dimension, |row, column| {
                e_matrix.at(row, column) - a_matrix.at(row, column) / a_norm_one
            }));
        let c_vector: Vec<Matrix::Scalar> =
            Vec::from_column(&ColumnFunc::new(dimension, |i| vector.at(i) / a_norm_one));

        return run_method(b_matrix, NormEnum::Two, c_vector, accuracy);
    }

    // multiply by A^T if none of the above work
    let new_a_matrix: DenseRowMatrix<Matrix::Scalar> = a_matrix.symmetric_square();
    let new_vector: Vec<Matrix::Scalar> = apply(a_matrix, vector);
    return internal_iterative_solve(&new_a_matrix, &new_vector, accuracy, true);
}

enum NormEnum<Scalar> {
    One { value: Scalar },
    Two,
    Infty { value: Scalar },
}

fn run_method<Matrix, ColumnIn, ColumnOut>(
    b_matrix: Matrix,
    b_norm: NormEnum<Matrix::Scalar>,
    c_vector: ColumnIn,
    accuracy: Matrix::Scalar,
) -> ColumnOut
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnMut<Scalar = Matrix::Scalar> + ColumnFuncInitializer<Scalar = Matrix::Scalar>,
{
    ColumnOut::new_fill(0, Matrix::Scalar::zero())
}
