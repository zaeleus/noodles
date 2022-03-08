//! SAM record sequence and bases.

pub mod base;

pub use self::base::Base;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
    str::FromStr,
};

use noodles_core::Position;

use super::NULL_FIELD;

/// A SAM record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence(Vec<Base>);

impl Sequence {
    /// Returns a reference to the base at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let sequence: Sequence = "ATCG".parse()?;
    ///
    /// let i = Position::try_from(2)?;
    /// assert_eq!(sequence.get(i), Some(&Base::T));
    ///
    /// let i = Position::try_from(8)?;
    /// assert!(sequence.get(i).is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get(&self, index: Position) -> Option<&Base> {
        let i = usize::from(index) - 1;
        self.0.get(i)
    }

    /// Returns a mutable reference to the base at the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let mut sequence: Sequence = "ATCG".parse()?;
    ///
    /// let i = Position::try_from(2)?;
    /// if let Some(base) = sequence.get_mut(i) {
    ///     *base = Base::N;
    /// }
    ///
    /// assert_eq!(sequence.get(i), Some(&Base::N));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get_mut(&mut self, index: Position) -> Option<&mut Base> {
        let i = usize::from(index) - 1;
        self.0.get_mut(i)
    }
}

impl Deref for Sequence {
    type Target = Vec<Base>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Sequence {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<Base>> for Sequence {
    fn from(bases: Vec<Base>) -> Self {
        Self(bases)
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            write!(f, "{}", NULL_FIELD)
        } else {
            for base in &self.0 {
                write!(f, "{}", base)?;
            }

            Ok(())
        }
    }
}

/// An error returned when a raw SAM record sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The raw sequence has an invalid base.
    InvalidBase(base::TryFromCharError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidBase(e) => write!(f, "invalid base: {}", e),
        }
    }
}

impl FromStr for Sequence {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            NULL_FIELD => Ok(Self::default()),
            _ => Self::try_from(s.as_bytes().to_vec()),
        }
    }
}

impl TryFrom<Vec<u8>> for Sequence {
    type Error = ParseError;

    fn try_from(buf: Vec<u8>) -> Result<Self, Self::Error> {
        buf.into_iter()
            .map(|b| b.to_ascii_uppercase())
            .map(Base::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map(Self::from)
            .map_err(ParseError::InvalidBase)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let sequence = Sequence::from(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!(sequence.to_string(), "ATCG");
    }

    #[test]
    fn test_from_str() {
        let expected = Sequence(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!("ATCG".parse::<Sequence>(), Ok(expected));

        let expected = Sequence(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!("atcg".parse::<Sequence>(), Ok(expected));

        let expected = Sequence(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!("aTcG".parse::<Sequence>(), Ok(expected));

        assert_eq!("*".parse::<Sequence>(), Ok(Sequence::default()));

        assert_eq!("".parse::<Sequence>(), Err(ParseError::Empty));
    }
}
