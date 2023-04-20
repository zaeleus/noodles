use std::{
    borrow::Borrow,
    error, fmt,
    hash::{Hash, Hasher},
    marker::PhantomData,
    str::FromStr,
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

/// An error returned when a raw VCF header record value map other tag fails to a parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl<S> FromStr for Other<S>
where
    S: Standard,
{
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use super::Tag;

        match Tag::from(s.to_string()) {
            Tag::Standard(_) => Err(ParseError::Invalid),
            Tag::Other(tag) => Ok(tag),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        use crate::header::record::value::map::tag::Identity;

        assert_eq!(
            "NOODLES".parse(),
            Ok(Other(String::from("NOODLES"), PhantomData::<Identity>))
        );
        assert_eq!("ID".parse::<Other<Identity>>(), Err(ParseError::Invalid));
    }
}
