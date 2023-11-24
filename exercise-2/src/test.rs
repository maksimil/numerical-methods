use crate::{
    basic::Index,
    matrix::{
        column::ColumnFunc,
        dense::DenseRowMatrix,
        norms::NormedColumn,
        traits::{MatrixFuncInitializer, MatrixRef},
    },
};

pub type Scalar = f64;

#[derive(Debug, Clone)]
pub struct TestCase {
    pub name: String,
    pub matrix: DenseRowMatrix<Scalar>,
    pub vector: Vec<Scalar>,
    pub answer: Vec<Scalar>,
}

#[derive(Debug, Clone)]
pub struct TestResult {
    pub answer: Vec<Scalar>,
    pub norm_one: Scalar,
}

pub fn on_case(
    case: &TestCase,
    implementation: impl Fn(&DenseRowMatrix<Scalar>, &Vec<Scalar>) -> Vec<Scalar>,
) -> TestResult {
    let dimension = case.matrix.dimension();
    let impl_answer = implementation(&case.matrix, &case.vector);
    let norm_one = ColumnFunc::new(dimension, |i| impl_answer[i] - case.answer[i]).norm_one();
    TestResult {
        answer: impl_answer,
        norm_one,
    }
}

pub fn create_fifth_case(epsilon: Scalar, dimension: Index) -> TestCase {
    let alpha = 6.0 * epsilon;

    let mut vector = vec![-1.0; dimension];
    vector[dimension - 1] = 1.0;

    let mut answer = vec![0.0; dimension];
    answer[dimension - 1] = 1.0 / (1.0 + alpha);

    let matrix = DenseRowMatrix::new_func(dimension, |i, j| {
        if i == j {
            1.0 + alpha
        } else if i > j {
            alpha
        } else {
            -1.0 - alpha
        }
    });

    TestCase {
        name: format!("5;{};{:e}", dimension, epsilon),
        matrix,
        answer,
        vector,
    }
}

pub fn create_fifth_cases(epsilons: Vec<Scalar>, dimensions: Vec<Index>) -> Vec<TestCase> {
    dimensions
        .into_iter()
        .map(|dimension| {
            epsilons
                .clone()
                .into_iter()
                .map(move |epsilon| create_fifth_case(epsilon, dimension))
        })
        .flatten()
        .collect()
}

pub fn create_static_test_cases() -> Vec<TestCase> {
    vec![
        TestCase {
            name: String::from("0"),
            matrix: DenseRowMatrix::new(3, vec![0.0, 2.0, 3.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0]),
            answer: vec![1.0, 2.0, 3.0],
            vector: vec![13.0, 17.0, 32.0],
        },
        TestCase {
            name: String::from("1"),
            matrix: DenseRowMatrix::new(3, vec![8.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 12.0]),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![1.0, 1.0, 1.0],
        },
        TestCase {
            name: String::from("2"),
            matrix: DenseRowMatrix::new(3, vec![-8.0, 1.0, 1.0, 1.0, -10.0, 1.0, 1.0, 1.0, -12.0]),
            vector: vec![-10.0, -12.0, -14.0],
            answer: vec![375.0 / 232.0, 349.0 / 232.0, 331.0 / 232.0],
        },
        TestCase {
            name: String::from("3"),
            matrix: DenseRowMatrix::new(
                3,
                vec![-8.0, 9.0, 10.0, 11.0, -10.0, 7.0, 10.0, 11.0, -12.0],
            ),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![444.0 / 307.0, 358.0 / 307.0, 340.0 / 307.0],
        },
        TestCase {
            name: String::from("4"),
            matrix: DenseRowMatrix::new(3, vec![8.0, 9.0, 10.0, 11.0, 10.0, 7.0, 10.0, 11.0, 12.0]),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![16.0, -22.0, 8.0],
        },
    ]
}
