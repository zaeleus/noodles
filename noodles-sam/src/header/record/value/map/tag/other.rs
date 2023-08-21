use std::{
    borrow::Borrow,
    fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
};

use super::LENGTH;

#[derive(Clone, Copy, Debug)]
pub struct Other<S>(pub(super) [u8; LENGTH], pub(super) PhantomData<S>);

impl<S> Borrow<[u8; LENGTH]> for Other<S> {
    fn borrow(&self) -> &[u8; LENGTH] {
        &self.0
    }
}

impl<S> Hash for Other<S> {
    fn hash<H>(&self, state: &mut H)
    where
        H: Hasher,
    {
        self.0.hash(state);
    }
}

impl<S> PartialEq for Other<S> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<S> Eq for Other<S> {}

impl<S> fmt::Display for Other<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        char::from(self.0[0]).fmt(f)?;
        char::from(self.0[1]).fmt(f)?;
        Ok(())
    }
}
