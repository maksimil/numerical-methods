use crate::{
    basic::Index,
    matrix::{
        dense::DenseRowMatrix,
        traits::{MatrixFuncInitializer, MatrixRef},
    },
};

type Scalar = f64;

#[derive(Debug, Clone)]
struct TestCase {
    name: String,
    matrix: DenseRowMatrix<Scalar>,
    vector: Vec<Scalar>,
    answer: Vec<Scalar>,
}

#[derive(Debug, Clone)]
struct TestResult {
    answer: Vec<Scalar>,
    norm_one: Scalar,
}

fn on_case(
    case: &TestCase,
    implementation: impl Fn(&DenseRowMatrix<Scalar>, &Vec<Scalar>) -> Vec<Scalar>,
) -> TestResult {
    let dimension = case.matrix.dimension();
    let impl_answer = implementation(&case.matrix, &case.vector);
    let norm_one = (0..dimension)
        .map(|i| (impl_answer[i] - case.answer[i]).abs())
        .sum();
    TestResult {
        answer: impl_answer,
        norm_one,
    }
}

fn create_fifth_case(epsilon: Scalar, dimension: Index) -> TestCase {
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
        name: format!("Test 5 with epsilon={:e}, n={}", epsilon, dimension),
        matrix,
        answer,
        vector,
    }
}

fn create_test_cases() -> Vec<TestCase> {
    let fifth_cases = vec![1e-3, 1e-6, 1e-9, 1e-12]
        .into_iter()
        .map(|epsilon| {
            vec![4, 5, 8, 16, 32]
                .into_iter()
                .map(move |dimension| create_fifth_case(epsilon, dimension))
        })
        .flatten();
    vec![
        TestCase {
            name: String::from("Test 0"),
            matrix: DenseRowMatrix::new(3, vec![0.0, 2.0, 3.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0]),
            answer: vec![1.0, 2.0, 3.0],
            vector: vec![13.0, 17.0, 32.0],
        },
        TestCase {
            name: String::from("Test 1"),
            matrix: DenseRowMatrix::new(3, vec![8.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 12.0]),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![1.0, 1.0, 1.0],
        },
        TestCase {
            name: String::from("Test 2"),
            matrix: DenseRowMatrix::new(3, vec![-8.0, 1.0, 1.0, 1.0, -10.0, 1.0, 1.0, 1.0, -12.0]),
            vector: vec![-10.0, -12.0, -14.0],
            answer: vec![375.0 / 232.0, 349.0 / 232.0, 331.0 / 232.0],
        },
        TestCase {
            name: String::from("Test 3"),
            matrix: DenseRowMatrix::new(
                3,
                vec![-8.0, 9.0, 10.0, 11.0, -10.0, 7.0, 10.0, 11.0, -12.0],
            ),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![444.0 / 307.0, 358.0 / 307.0, 340.0 / 307.0],
        },
        TestCase {
            name: String::from("Test 4"),
            matrix: DenseRowMatrix::new(3, vec![8.0, 9.0, 10.0, 11.0, 10.0, 7.0, 10.0, 11.0, 12.0]),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![16.0, -22.0, 8.0],
        },
    ]
    .into_iter()
    .chain(fifth_cases)
    .collect()
}

pub fn test_cases(
    implementation: impl Fn(&DenseRowMatrix<Scalar>, &Vec<Scalar>) -> Vec<Scalar> + Clone,
) {
    for case in create_test_cases() {
        let result = on_case(&case, implementation.clone());
        println!("name={}, norm_inf={:e}", case.name, result.norm_one);
    }
}
