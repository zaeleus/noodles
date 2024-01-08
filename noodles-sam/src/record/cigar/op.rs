//! SAM record CIGAR operation and kind.

pub mod kind;

use std::fmt;

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
    pub const fn new(kind: Kind, len: usize) -> Self {
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
}
