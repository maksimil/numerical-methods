use crate::{lu_decomposition::lu_decomposition::LUDecomposition, test::test_cases};

mod basic;
mod lu_decomposition;
mod matrix;
mod representation;
mod test;

fn main() {
    test_cases(|matrix, vector| {
        let decomposition = LUDecomposition::calculate(matrix.clone());
        let mut vector_mut = vector.clone();
        decomposition.solve(&mut vector_mut);
        vector_mut
    });
}
