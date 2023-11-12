use crate::basic::Index;

pub trait MatrixRef<Scalar> {
    fn dimension(&self) -> Index;
    fn at(&self, row: Index, column: Index) -> &Scalar;
}

pub trait MatrixMutRef<Scalar>: MatrixRef<Scalar> {
    fn at_mut(&mut self, row: Index, column: Index) -> &mut Scalar;
}

pub trait MatrixFuncInitializer<Scalar> {
    fn new_func(dimension: Index, fill: impl Fn(Index, Index) -> Scalar) -> Self;
}

pub trait MatrixFillInitializer<Scalar> {
    fn new_fill(dimension: Index, fill: Scalar) -> Self;
}

impl<Scalar: Clone, T: MatrixFuncInitializer<Scalar>> MatrixFillInitializer<Scalar> for T {
    fn new_fill(dimension: Index, fill: Scalar) -> Self {
        T::new_func(dimension, |_, _| fill.clone())
    }
}
