//! Alignment record sequence and bases.

pub mod base;

pub use self::base::Base;

use std::{
    error, fmt,
    ops::{Index, IndexMut},
    str::FromStr,
};

use noodles_core::position::SequenceIndex;

/// An alignment record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence(Vec<Base>);

impl Sequence {
    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Sequence;
    /// let sequence = Sequence::default();
    /// assert_eq!(sequence.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether the sequence is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Sequence;
    /// let sequence = Sequence::default();
    /// assert!(sequence.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Removes all bases from the sequence.
    ///
    /// This does not affect the capacity of this sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let mut sequence = Sequence::from(vec![Base::N]);
    /// assert!(!sequence.is_empty());
    ///
    /// sequence.clear();
    /// assert!(sequence.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }

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
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: SequenceIndex<Base>,
    {
        index.get(self.0.as_ref())
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
    pub fn get_mut<I>(&mut self, index: I) -> Option<&mut I::Output>
    where
        I: SequenceIndex<Base>,
    {
        index.get_mut(self.0.as_mut())
    }

    /// Appends a base to the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{sequence::Base, Sequence};
    ///
    /// let mut sequence = Sequence::default();
    /// sequence.push(Base::N);
    ///
    /// let expected = Sequence::from(vec![Base::N]);
    ///
    /// assert_eq!(sequence, expected);
    /// ````
    pub fn push(&mut self, base: Base) {
        self.0.push(base);
    }
}

impl AsRef<[Base]> for Sequence {
    fn as_ref(&self) -> &[Base] {
        &self.0
    }
}

impl AsMut<Vec<Base>> for Sequence {
    fn as_mut(&mut self) -> &mut Vec<Base> {
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
        for base in self.as_ref() {
            write!(f, "{}", base)?;
        }

        Ok(())
    }
}

/// An error returned when a raw alignment record sequence fails to parse.
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
        if s.is_empty() {
            Err(ParseError::Empty)
        } else {
            Self::try_from(s.as_bytes().to_vec())
        }
    }
}

impl<I> Index<I> for Sequence
where
    I: SequenceIndex<Base>,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        index.index(&self.0)
    }
}

impl<I> IndexMut<I> for Sequence
where
    I: SequenceIndex<Base>,
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        index.index_mut(&mut self.0)
    }
}

impl TryFrom<Vec<u8>> for Sequence {
    type Error = ParseError;

    fn try_from(buf: Vec<u8>) -> Result<Self, Self::Error> {
        buf.into_iter()
            .map(Base::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map(Self::from)
            .map_err(ParseError::InvalidBase)
    }
}

impl From<Sequence> for Vec<Base> {
    fn from(sequence: Sequence) -> Self {
        sequence.0
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

        assert_eq!("".parse::<Sequence>(), Err(ParseError::Empty));

        assert!(matches!(
            "*".parse::<Sequence>(),
            Err(ParseError::InvalidBase(_))
        ));
    }
}
