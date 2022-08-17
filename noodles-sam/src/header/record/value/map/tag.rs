use std::{
    borrow::Borrow,
    fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
    str::FromStr,
};

pub(crate) const LENGTH: usize = 2;

pub trait Standard: TryFrom<[u8; LENGTH]> {}

#[derive(Clone, Debug)]
pub struct Other<S>([u8; LENGTH], PhantomData<S>);

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

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Tag<S> {
    Standard(S),
    Other(Other<S>),
}

impl<S> From<[u8; LENGTH]> for Tag<S>
where
    S: Standard,
{
    fn from(s: [u8; LENGTH]) -> Self {
        match S::try_from(s) {
            Ok(tag) => Self::Standard(tag),
            Err(_) => Self::Other(Other(s, PhantomData)),
        }
    }
}

impl<S> FromStr for Tag<S>
where
    S: Standard,
{
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let b: [u8; 2] = s.as_bytes().try_into().map_err(|_| ())?;
        Ok(Self::from(b))
    }
}
