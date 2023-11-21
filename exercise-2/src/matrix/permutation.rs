use std::borrow::Borrow;

use crate::{basic::Index, representation::Representation};

pub trait Permutation {
    fn permute(&self, index: Index) -> Index;
}

impl Permutation for Vec<Index> {
    fn permute(&self, index: Index) -> Index {
        self[index]
    }
}

impl<Stored, Impl> Permutation for Representation<Stored, Impl>
where
    Stored: Borrow<Impl>,
    Impl: Permutation,
{
    fn permute(&self, index: Index) -> Index {
        self.represent().permute(index)
    }
}
