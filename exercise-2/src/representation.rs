use std::{
    borrow::{Borrow, BorrowMut},
    marker::PhantomData,
};

// An owned version of Borrow and BorrowMut
#[derive(Clone, Debug)]
pub struct Representation<Stored, Impl> {
    stored: Stored,
    _phantom: PhantomData<Impl>,
}

impl<Stored, Impl> From<Stored> for Representation<Stored, Impl> {
    fn from(stored: Stored) -> Self {
        Representation {
            stored,
            _phantom: PhantomData {},
        }
    }
}

impl<Stored, Impl> Representation<Stored, Impl> {
    pub fn represent(&self) -> &Impl
    where
        Stored: Borrow<Impl>,
    {
        self.stored.borrow()
    }

    pub fn represent_mut(&mut self) -> &mut Impl
    where
        Stored: BorrowMut<Impl>,
    {
        self.stored.borrow_mut()
    }
}

pub type RefRepresentation<'a, T> = Representation<&'a T, T>;
pub type MutRefRepresentation<'a, T> = Representation<&'a mut T, T>;

pub fn repr_ref<'a, T>(v: &'a T) -> RefRepresentation<'a, T> {
    RefRepresentation::from(v)
}

pub fn repr_mut<'a, T>(v: &'a mut T) -> MutRefRepresentation<'a, T> {
    MutRefRepresentation::from(v)
}
