//! SAM record CIGAR operation and kind.

pub mod kind;

use std::{error, fmt, num, str::FromStr};

pub use self::kind::Kind;

/// A SAM record CIGAR operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Op {
    kind: Kind,
    len: usize,
}

impl Op {
    /// Creates a CIGAR operation.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::cigar::{op::Kind, Op};
    ///
    /// let op = Op::new(Kind::Match, 13);
    ///
    /// assert_eq!(op.kind(), Kind::Match);
    /// assert_eq!(op.len(), 13);
    /// ```
    pub fn new(kind: Kind, len: usize) -> Self {
        Self { kind, len }
    }

    /// Returns the kind of the operation.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::cigar::{op::Kind, Op};
    /// let op = Op::new(Kind::Match, 13);
    /// assert_eq!(op.kind(), Kind::Match);
    /// ```
    pub fn kind(self) -> Kind {
        self.kind
    }

    /// Returns the length of the operation.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::cigar::{op::Kind, Op};
    /// let op = Op::new(Kind::Match, 13);
    /// assert_eq!(op.len(), 13);
    /// ```
    pub fn len(self) -> usize {
        self.len
    }

    /// Returns whether the operation is a no-op.
    ///
    /// That is, whether the operation has a length of 0.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::cigar::{op::Kind, Op};
    ///
    /// let op = Op::new(Kind::Match, 0);
    /// assert!(op.is_empty());
    ///
    /// let op = Op::new(Kind::Match, 13);
    /// assert!(!op.is_empty());
    /// ```
    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

impl fmt::Display for Op {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.kind())
    }
}

/// An error returned when a raw CIGAR operation fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The raw operation has an invalid length.
    InvalidLength(num::ParseIntError),
    /// The raw operation has an invalid kind.
    InvalidKind(kind::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidLength(e) => write!(f, "invalid length: {}", e),
            Self::InvalidKind(e) => write!(f, "invalid kind: {}", e),
        }
    }
}

impl FromStr for Op {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let (raw_len, raw_kind) = s.split_at(s.len() - 1);

        let len = raw_len.parse().map_err(ParseError::InvalidLength)?;
        let kind = raw_kind.parse().map_err(ParseError::InvalidKind)?;

        Ok(Self::new(kind, len))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let op = Op::new(Kind::Match, 5);
        assert_eq!(op.to_string(), "5M");

        let op = Op::new(Kind::SoftClip, 13);
        assert_eq!(op.to_string(), "13S");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("1M".parse(), Ok(Op::new(Kind::Match, 1)));
        assert_eq!("13N".parse(), Ok(Op::new(Kind::Skip, 13)));
        assert_eq!("144S".parse(), Ok(Op::new(Kind::SoftClip, 144)));

        assert_eq!("".parse::<Op>(), Err(ParseError::Empty));

        assert!(matches!(
            "Z".parse::<Op>(),
            Err(ParseError::InvalidLength(_))
        ));

        assert!(matches!(
            "21".parse::<Op>(),
            Err(ParseError::InvalidKind(_))
        ));

        Ok(())
    }
}
