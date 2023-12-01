use crate::{
    basic::{Index, Numerical, OtherNumericalOps},
    matrix::{
        column::{
            apply, apply_at, dot, ColumnFunc, ColumnFuncInitializer, ColumnMut, ColumnOf, ColumnRef,
        },
        dense::DenseRowMatrix,
        norms::{NormedColumn, NormedMatrix},
        traits::{MatrixFunc, MatrixFuncInitializer, MatrixRef},
        transpose::{MatrixSymmetricSquare, MatrixTranspose},
    },
    representation::repr_ref,
};

pub fn simple_iterative_solve<Matrix, ColumnIn, ColumnOut>(
    matrix: &Matrix,
    vector: &ColumnIn,
    accuracy: Matrix::Scalar,
) -> (ColumnOut, Index)
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnMut<Scalar = Matrix::Scalar> + ColumnFuncInitializer<Scalar = Matrix::Scalar>,
{
    internal_simple_iterative_solve(matrix, vector, accuracy, false)
}

fn internal_simple_iterative_solve<Matrix, ColumnIn, ColumnOut>(
    a_matrix: &Matrix,
    vector: &ColumnIn,
    accuracy: Matrix::Scalar,
    can_abuse_norm: bool,
) -> (ColumnOut, Index)
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

        // println!("norm one: |A|={:?} |B|={:?}", a_norm_one, q_one);

        if q_one < Matrix::Scalar::one() {
            let b_matrix: DenseRowMatrix<Matrix::Scalar> = DenseRowMatrix::from_matrix(&b_matrix);
            let c_vector = ColumnOut::new_func(dimension, |i| vector.at(i) / a_norm_one);
            return run_method(
                a_matrix,
                vector,
                &b_matrix,
                NormEnum::One { value: q_one },
                &c_vector,
                accuracy,
            );
        }
    }

    // try p=inf
    {
        let a_norm_inf = a_matrix.norm_inf();
        let b_matrix = MatrixFunc::new(dimension, |row, column| {
            e_matrix.at(row, column) - a_matrix.at(row, column) / a_norm_inf
        });
        let q_inf = b_matrix.norm_inf();

        // println!("norm inf: |A|={:?} |B|={:?}", a_norm_inf, q_inf);

        if q_inf < Matrix::Scalar::one() {
            let b_matrix: DenseRowMatrix<Matrix::Scalar> = DenseRowMatrix::from_matrix(&b_matrix);
            let c_vector = ColumnOut::new_func(dimension, |i| vector.at(i) / a_norm_inf);
            return run_method(
                a_matrix,
                vector,
                &b_matrix,
                NormEnum::Infty { value: q_inf },
                &c_vector,
                accuracy,
            );
        }
    }

    // abuse the p=2 norm if can
    if can_abuse_norm {
        let a_norm_one = a_matrix.norm_one();
        let b_matrix: DenseRowMatrix<Matrix::Scalar> =
            DenseRowMatrix::from_matrix(&MatrixFunc::new(dimension, |row, column| {
                e_matrix.at(row, column) - a_matrix.at(row, column) / a_norm_one
            }));
        let c_vector = ColumnOut::new_func(dimension, |i| vector.at(i) / a_norm_one);

        // println!("Abusing the second norm");

        return run_method(
            a_matrix,
            vector,
            &b_matrix,
            NormEnum::Two,
            &c_vector,
            accuracy,
        );
    }

    // multiply by A^T if none of the above work
    let new_a_matrix: DenseRowMatrix<Matrix::Scalar> = a_matrix.symmetric_square();
    let new_vector: Vec<Matrix::Scalar> =
        apply(&MatrixTranspose::from(repr_ref::<Matrix>(a_matrix)), vector);
    return internal_simple_iterative_solve(&new_a_matrix, &new_vector, accuracy, true);
}

pub fn zeidel_iterative_solve<Matrix, ColumnIn, ColumnOut>(
    matrix: &Matrix,
    vector: &ColumnIn,
    accuracy: Matrix::Scalar,
) -> (ColumnOut, Index)
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnMut<Scalar = Matrix::Scalar> + ColumnFuncInitializer<Scalar = Matrix::Scalar>,
{
    internal_zeidel_iterative_solve(matrix, vector, accuracy, false)
}

fn internal_zeidel_iterative_solve<Matrix, ColumnIn, ColumnOut>(
    a_matrix: &Matrix,
    vector: &ColumnIn,
    accuracy: Matrix::Scalar,
    can_abuse_norm: bool,
) -> (ColumnOut, Index)
where
    Matrix: MatrixRef,
    Matrix::Scalar: Numerical,
    ColumnIn: ColumnRef<Scalar = Matrix::Scalar>,
    ColumnOut: ColumnMut<Scalar = Matrix::Scalar> + ColumnFuncInitializer<Scalar = Matrix::Scalar>,
{
    let dimension = a_matrix.dimension();

    // checking with diagonal
    if (0..dimension)
        .map(|i| a_matrix.at(i, i) != Matrix::Scalar::zero())
        .min()
        .unwrap()
    {
        // construct the b matrix
        let b_matrix = DenseRowMatrix::new_func(dimension, |row, column| {
            if row == column {
                Matrix::Scalar::zero()
            } else {
                -a_matrix.at(row, column) / a_matrix.at(row, row)
            }
        });

        let c_vector = ColumnOut::new_func(dimension, |i| vector.at(i) / a_matrix.at(i, i));

        // testing p=1
        let q_one = b_matrix.norm_one();
        if q_one < Matrix::Scalar::one() {
            return run_zeidel_method(
                a_matrix,
                vector,
                &b_matrix,
                NormEnum::One { value: q_one },
                &c_vector,
                accuracy,
            );
        }

        // testing p=inf
        let q_inf = b_matrix.norm_inf();
        if q_inf < Matrix::Scalar::one() {
            return run_zeidel_method(
                a_matrix,
                vector,
                &b_matrix,
                NormEnum::Infty { value: q_inf },
                &c_vector,
                accuracy,
            );
        }

        if can_abuse_norm {
            return run_zeidel_method(
                a_matrix,
                vector,
                &b_matrix,
                NormEnum::Two,
                &c_vector,
                accuracy,
            );
        }
    }

    // multiply by A^T if none of the above work
    let new_a_matrix: DenseRowMatrix<Matrix::Scalar> = a_matrix.symmetric_square();
    let new_vector: Vec<Matrix::Scalar> =
        apply(&MatrixTranspose::from(repr_ref::<Matrix>(a_matrix)), vector);
    return internal_zeidel_iterative_solve(&new_a_matrix, &new_vector, accuracy, true);
}

#[derive(Clone)]
enum NormEnum<Scalar> {
    One { value: Scalar },
    Two,
    Infty { value: Scalar },
}

fn run_method<AMatrix, BMatrix, ColumnIn, ColumnOut>(
    a_matrix: &AMatrix,
    vector: &ColumnIn,
    b_matrix: &BMatrix,
    b_norm: NormEnum<BMatrix::Scalar>,
    c_vector: &ColumnOut,
    accuracy: BMatrix::Scalar,
) -> (ColumnOut, Index)
where
    BMatrix: MatrixRef,
    BMatrix::Scalar: Numerical,
    AMatrix: MatrixRef<Scalar = BMatrix::Scalar>,
    ColumnIn: ColumnRef<Scalar = BMatrix::Scalar>,
    ColumnOut:
        ColumnMut<Scalar = BMatrix::Scalar> + ColumnFuncInitializer<Scalar = BMatrix::Scalar>,
{
    let mut steps = 0;
    // println!(
    //     "Running B={:#?} x0={:#?}",
    //     DenseRowMatrix::<BMatrix::Scalar>::from_matrix(b_matrix),
    //     Vec::<BMatrix::Scalar>::from_column(c_vector)
    // );
    let dimension = b_matrix.dimension();
    let mut x_this = ColumnOut::new_func(dimension, |i| c_vector.at(i));
    let mut achieved_accuracy = accuracy;
    while achieved_accuracy >= accuracy {
        let mut x_next = ColumnOut::new_func(dimension, |i| {
            apply_at(b_matrix, &x_this, i) + c_vector.at(i)
        });
        // println!("x_k+1={:#?}", Vec::<BMatrix::Scalar>::from_column(&x_next));
        achieved_accuracy = calculate_accuracy(a_matrix, vector, b_norm.clone(), &x_this, &x_next);
        // println!("Acc={:?}", achieved_accuracy);
        std::mem::swap(&mut x_this, &mut x_next);
        steps += 1;
    }
    (x_this, steps)
}

fn run_zeidel_method<AMatrix, BMatrix, ColumnIn, ColumnOut>(
    a_matrix: &AMatrix,
    vector: &ColumnIn,
    b_matrix: &BMatrix,
    b_norm: NormEnum<BMatrix::Scalar>,
    c_vector: &ColumnOut,
    accuracy: BMatrix::Scalar,
) -> (ColumnOut, Index)
where
    BMatrix: MatrixRef,
    BMatrix::Scalar: Numerical,
    AMatrix: MatrixRef<Scalar = BMatrix::Scalar>,
    ColumnIn: ColumnRef<Scalar = BMatrix::Scalar>,
    ColumnOut:
        ColumnMut<Scalar = BMatrix::Scalar> + ColumnFuncInitializer<Scalar = BMatrix::Scalar>,
{
    let mut steps = 0;
    let dimension = b_matrix.dimension();
    let mut x_this = ColumnOut::new_func(dimension, |i| c_vector.at(i));
    let mut achieved_accuracy = accuracy;
    while achieved_accuracy >= accuracy {
        let mut x_next = ColumnOut::new_func(dimension, |i| x_this.at(i));
        for i in 0..dimension {
            *x_next.at_mut(i) = dot(
                &ColumnOf::new(MatrixTranspose::from(repr_ref::<BMatrix>(b_matrix)), i),
                &x_next,
            ) + c_vector.at(i);
        }
        // println!("x_k+1={:#?}", Vec::<BMatrix::Scalar>::from_column(&x_next));
        achieved_accuracy = calculate_accuracy(a_matrix, vector, b_norm.clone(), &x_this, &x_next);
        // println!("Acc={:?}", achieved_accuracy);
        std::mem::swap(&mut x_this, &mut x_next);
        steps += 1;
    }
    (x_this, steps)
}

fn calculate_accuracy<Matrix, Column, ColumnIn>(
    a_matrix: &Matrix,
    vector: &ColumnIn,
    b_norm: NormEnum<Column::Scalar>,
    x_this: &Column,
    x_next: &Column,
) -> Column::Scalar
where
    Column: ColumnRef,
    Column::Scalar: Numerical,
    Matrix: MatrixRef<Scalar = Column::Scalar>,
    ColumnIn: ColumnRef<Scalar = Column::Scalar>,
{
    match b_norm {
        NormEnum::One { value } => {
            value / (Column::Scalar::one() - value)
                * ColumnFunc::new(x_this.dimension(), |i| x_this.at(i) - x_next.at(i)).norm_one()
        }
        NormEnum::Infty { value } => {
            value / (Column::Scalar::one() - value)
                * ColumnFunc::new(x_this.dimension(), |i| x_this.at(i) - x_next.at(i)).norm_inf()
        }
        NormEnum::Two => ColumnFunc::new(x_next.dimension(), |i| {
            apply_at(a_matrix, x_next, i) - vector.at(i)
        })
        .norm_one(),
    }
}
