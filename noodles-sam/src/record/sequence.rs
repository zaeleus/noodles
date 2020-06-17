//! SAM record sequence and bases.

mod base;

pub use self::base::Base;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use super::NULL_FIELD;

/// A SAM record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence {
    bases: Vec<Base>,
}

impl Sequence {
    /// Creates a new SAM record sequence from a list of bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{sequence::Base, Sequence};
    /// let sequence = Sequence::new(vec![Base::A, Base::C, Base::G, Base::T]);
    /// ```
    pub fn new(bases: Vec<Base>) -> Self {
        Self { bases }
    }

    /// Returns the list of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let sequence = Sequence::new(vec![Base::A, Base::C, Base::G, Base::T]);
    ///
    /// let actual = sequence.bases();
    /// let expected = [Base::A, Base::C, Base::G, Base::T];
    /// assert_eq!(actual, expected);
    /// ```
    pub fn bases(&self) -> &[Base] {
        &self.bases
    }

    /// Returns whether the sequence is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let sequence = Sequence::default();
    /// assert!(sequence.is_empty());
    ///
    /// let sequence = Sequence::new(vec![Base::A, Base::C, Base::G, Base::T]);
    /// assert!(!sequence.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let sequence = Sequence::default();
    /// assert_eq!(sequence.len(), 0);
    ///
    /// let sequence = Sequence::new(vec![Base::A, Base::C, Base::G, Base::T]);
    /// assert_eq!(sequence.len(), 4);
    /// ```
    pub fn len(&self) -> usize {
        self.bases.len()
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.bases.is_empty() {
            write!(f, "{}", NULL_FIELD)
        } else {
            for base in &self.bases {
                write!(f, "{}", base)?;
            }

            Ok(())
        }
    }
}

/// An error returned when a raw SAM record sequence fails to parse.
#[derive(Debug)]
pub struct ParseError(base::TryFromCharError);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid sequence: {}", self.0)
    }
}

impl FromStr for Sequence {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == NULL_FIELD {
            return Ok(Self::default());
        }

        s.chars()
            .map(|c| c.to_ascii_uppercase())
            .map(Base::try_from)
            .collect::<Result<_, _>>()
            .map(Self::new)
            .map_err(ParseError)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bases() {
        let bases = [Base::A, Base::T, Base::C, Base::G];
        let sequence = Sequence::new(bases.to_vec());
        assert_eq!(sequence.bases(), &bases);
    }

    #[test]
    fn test_is_empty() {
        let sequence = Sequence::new(Vec::new());
        assert!(sequence.is_empty());
    }

    #[test]
    fn test_len() {
        let sequence = Sequence::new(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!(sequence.len(), 4);
    }

    #[test]
    fn test_default() {
        assert!(Sequence::default().is_empty());
    }

    #[test]
    fn test_fmt() {
        let sequence = Sequence::new(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!(sequence.to_string(), "ATCG");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let expected = [Base::A, Base::T, Base::C, Base::G];

        let sequence = "ATCG".parse::<Sequence>()?;
        assert_eq!(sequence.bases(), &expected);

        let sequence = "atcg".parse::<Sequence>()?;
        assert_eq!(sequence.bases(), &expected);

        let sequence = "aTcG".parse::<Sequence>()?;
        assert_eq!(sequence.bases(), &expected);

        Ok(())
    }
}
