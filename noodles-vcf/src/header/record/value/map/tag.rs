use std::{
    borrow::Borrow,
    fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
    str::FromStr,
};

pub(super) const ID: &str = "ID";
pub(super) const NUMBER: &str = "Number";
pub(super) const TYPE: &str = "Type";
pub(super) const DESCRIPTION: &str = "Description";
pub(super) const IDX: &str = "IDX";

pub trait Standard: AsRef<str> + FromStr {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Identity {
    Id,
}

impl Standard for Identity {}

impl AsRef<str> for Identity {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => ID,
        }
    }
}

impl FromStr for Identity {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ID => Ok(Self::Id),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TypedDescribedIndexed {
    Id,
    Number,
    Type,
    Description,
    Idx,
}

impl Standard for TypedDescribedIndexed {}

impl AsRef<str> for TypedDescribedIndexed {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => ID,
            Self::Number => NUMBER,
            Self::Type => TYPE,
            Self::Description => DESCRIPTION,
            Self::Idx => IDX,
        }
    }
}

impl FromStr for TypedDescribedIndexed {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ID => Ok(Self::Id),
            NUMBER => Ok(Self::Number),
            TYPE => Ok(Self::Type),
            DESCRIPTION => Ok(Self::Description),
            IDX => Ok(Self::Idx),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Described {
    Id,
    Description,
}

impl Standard for Described {}

impl AsRef<str> for Described {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => ID,
            Self::Description => DESCRIPTION,
        }
    }
}

impl FromStr for Described {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ID => Ok(Self::Id),
            DESCRIPTION => Ok(Self::Description),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DescribedIndexed {
    Id,
    Description,
    Idx,
}

impl Standard for DescribedIndexed {}

impl AsRef<str> for DescribedIndexed {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => ID,
            Self::Description => DESCRIPTION,
            Self::Idx => IDX,
        }
    }
}

impl FromStr for DescribedIndexed {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ID => Ok(Self::Id),
            DESCRIPTION => Ok(Self::Description),
            IDX => Ok(Self::Idx),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Other<S>(String, PhantomData<S>);

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

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Tag<S> {
    Standard(S),
    Other(Other<S>),
}

impl<S> AsRef<str> for Tag<S>
where
    S: Standard,
{
    fn as_ref(&self) -> &str {
        match self {
            Self::Standard(tag) => tag.as_ref(),
            Self::Other(tag) => tag.as_ref(),
        }
    }
}

impl<S> From<String> for Tag<S>
where
    S: Standard,
{
    fn from(s: String) -> Self {
        match s.parse::<S>() {
            Ok(tag) => Self::Standard(tag),
            Err(_) => Self::Other(Other(s, PhantomData)),
        }
    }
}
