use crate::matrix::{dense::DenseRowMatrix, traits::MatrixRef};

type Scalar = f64;

#[derive(Debug, Clone)]
struct TestCase {
    matrix: DenseRowMatrix<Scalar>,
    vector: Vec<Scalar>,
    answer: Vec<Scalar>,
}

#[derive(Debug, Clone)]
struct TestResult {
    answer: Vec<Scalar>,
    norm_inf: Scalar,
}

fn on_case(
    case: &TestCase,
    implementation: impl Fn(&DenseRowMatrix<Scalar>, &Vec<Scalar>) -> Vec<Scalar>,
) -> TestResult {
    let dimension = case.matrix.dimension();
    let impl_answer = implementation(&case.matrix, &case.vector);
    let norm_inf = (0..dimension)
        .map(|i| (impl_answer[i] - case.answer[i]).abs())
        .fold(0.0, |m, v| if v > m { v } else { m });
    TestResult {
        answer: impl_answer,
        norm_inf,
    }
}

fn create_test_cases() -> Vec<TestCase> {
    vec![
        TestCase {
            matrix: DenseRowMatrix::new(3, vec![0.0, 2.0, 3.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0]),
            answer: vec![1.0, 2.0, 3.0],
            vector: vec![13.0, 17.0, 32.0],
        },
        TestCase {
            matrix: DenseRowMatrix::new(3, vec![8.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 12.0]),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![1.0, 1.0, 1.0],
        },
        TestCase {
            matrix: DenseRowMatrix::new(3, vec![-8.0, 1.0, 1.0, 1.0, -10.0, 1.0, 1.0, 1.0, -12.0]),
            vector: vec![-10.0, -12.0, -14.0],
            answer: vec![0.0, 0.0, 0.0],
        },
        TestCase {
            matrix: DenseRowMatrix::new(
                3,
                vec![-8.0, 9.0, 10.0, 11.0, -10.0, 7.0, 10.0, 11.0, -12.0],
            ),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![0.0, 0.0, 0.0],
        },
        TestCase {
            matrix: DenseRowMatrix::new(3, vec![8.0, 9.0, 10.0, 11.0, 10.0, 7.0, 10.0, 11.0, 12.0]),
            vector: vec![10.0, 12.0, 14.0],
            answer: vec![0.0, 0.0, 0.0],
        },
    ]
}

pub fn test_cases(
    implementation: impl Fn(&DenseRowMatrix<Scalar>, &Vec<Scalar>) -> Vec<Scalar> + Clone,
) {
    for case in create_test_cases() {
        let result = on_case(&case, implementation.clone());
        println!("Result = {:?}", result);
    }
}
