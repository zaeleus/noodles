use std::{
    borrow::Borrow,
    fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
};

use super::Standard;

#[derive(Clone, Debug)]
pub struct Other<S>(pub(super) String, pub(super) PhantomData<S>);

impl<S> AsRef<str> for Other<S>
where
    S: Standard,
{
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl<S> Borrow<str> for Other<S> {
    fn borrow(&self) -> &str {
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
        self.0.fmt(f)
    }
}
