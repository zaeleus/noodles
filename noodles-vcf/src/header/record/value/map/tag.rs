use std::{
    borrow::Borrow,
    fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
    str::FromStr,
};

pub trait Standard: FromStr {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Identity {
    Id,
}

impl Standard for Identity {}

impl FromStr for Identity {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
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

impl FromStr for TypedDescribedIndexed {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "Number" => Ok(Self::Number),
            "Type" => Ok(Self::Type),
            "Description" => Ok(Self::Description),
            "IDX" => Ok(Self::Idx),
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

impl FromStr for Described {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "Description" => Ok(Self::Description),
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

impl FromStr for DescribedIndexed {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "Description" => Ok(Self::Description),
            "IDX" => Ok(Self::Idx),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Other<S>(String, PhantomData<S>);

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
