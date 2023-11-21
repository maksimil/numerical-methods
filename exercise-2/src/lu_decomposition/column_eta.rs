use crate::basic::{Index, INDEX_NOT_FOUND};

#[derive(Debug, Clone)]
pub struct ColumnEtaMatrix<Scalar> {
    pub column: Index,
    pub eta_vector: Vec<Scalar>,
}

impl<Scalar> ColumnEtaMatrix<Scalar> {
    pub fn new(dimension: Index) -> ColumnEtaMatrix<Scalar>
    where
        Scalar: Default + Clone,
    {
        ColumnEtaMatrix {
            column: INDEX_NOT_FOUND,
            eta_vector: vec![Scalar::default(); dimension],
        }
    }
}
