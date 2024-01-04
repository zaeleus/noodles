//! SAM CIGAR and operations.

pub mod op;

use std::{error, fmt, io, ops::Deref, str::FromStr};

use self::op::Kind;
pub use self::op::Op;

/// A SAM record CIGAR.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Cigar(Vec<Op>);

impl Cigar {
    /// Removes all operations from the CIGAR.
    ///
    /// This does not affect its capacity.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{cigar::{op::Kind, Op}, Cigar};
    ///
    /// let mut cigar: Cigar = [Op::new(Kind::Match, 5)].into_iter().collect();
    /// assert!(!cigar.is_empty());
    ///
    /// cigar.clear();
    /// assert!(cigar.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Calculates the alignment span over the reference sequence.
    ///
    /// This sums the lengths of the CIGAR operations that consume the reference sequence, i.e.,
    /// alignment matches (`M`), deletions from the reference (`D`), skipped reference regions
    /// (`S`), sequence matches (`=`), and sequence mismatches (`X`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{cigar::{op::Kind, Op}, Cigar};
    ///
    /// let cigar: Cigar = [
    ///     Op::new(Kind::Match, 36),
    ///     Op::new(Kind::Deletion, 4),
    ///     Op::new(Kind::SoftClip, 8),
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// assert_eq!(cigar.alignment_span(), 40);
    /// ```
    pub fn alignment_span(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| op.kind().consumes_reference().then_some(op.len()))
            .sum()
    }

    /// Calculates the read length.
    ///
    /// This sums the lengths of the CIGAR operations that consume the read, i.e., alignment
    /// matches (`M`), insertions to the reference (`I`), soft clips (`S`), sequence matches (`=`),
    /// and sequence mismatches (`X`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{cigar::{op::Kind, Op}, Cigar};
    ///
    /// let cigar: Cigar = [
    ///     Op::new(Kind::Match, 36),
    ///     Op::new(Kind::Deletion, 4),
    ///     Op::new(Kind::SoftClip, 8),
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// assert_eq!(cigar.read_length(), 44);
    /// ```
    pub fn read_length(&self) -> usize {
        self.0
            .iter()
            .filter_map(|op| op.kind().consumes_read().then_some(op.len()))
            .sum()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Kind, usize)>> + '_> {
        Box::new(self.0.iter().map(|op| Ok((op.kind(), op.len()))))
    }
}

impl crate::alignment::record::Cigar for Cigar {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Kind, usize)>> + '_> {
        self.iter()
    }
}

impl crate::alignment::record::Cigar for &Cigar {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Kind, usize)>> + '_> {
        Cigar::iter(self)
    }
}

impl Deref for Cigar {
    type Target = [Op];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsMut<Vec<Op>> for Cigar {
    fn as_mut(&mut self) -> &mut Vec<Op> {
        &mut self.0
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for op in self.0.iter() {
            write!(f, "{op}")?;
        }

        Ok(())
    }
}

impl Extend<Op> for Cigar {
    fn extend<T: IntoIterator<Item = Op>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<Op> for Cigar {
    fn from_iter<T: IntoIterator<Item = Op>>(iter: T) -> Self {
        let mut cigar = Cigar::default();
        cigar.extend(iter);
        cigar
    }
}

impl TryFrom<Vec<Op>> for Cigar {
    type Error = ParseError;

    fn try_from(ops: Vec<Op>) -> Result<Self, Self::Error> {
        Ok(Self(ops))
    }
}

/// An error returned when a raw CIGAR string fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The CIGAR string has an invalid operation.
    InvalidOp(op::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidOp(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidOp(_) => f.write_str("invalid op"),
        }
    }
}

impl FromStr for Cigar {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut ops = Vec::new();

        let matches = s.match_indices(|c: char| !c.is_ascii_digit());
        let mut start = 0;

        for (end, raw_kind) in matches {
            let op = s[start..=end].parse().map_err(ParseError::InvalidOp)?;
            ops.push(op);
            start = end + raw_kind.len();
        }

        if start == s.len() {
            Self::try_from(ops)
        } else {
            Err(ParseError::Invalid)
        }
    }
}

impl From<Cigar> for Vec<Op> {
    fn from(cigar: Cigar) -> Self {
        cigar.0
    }
}

#[cfg(test)]
mod tests {
    use super::{op::Kind, *};

    #[test]
    fn test_is_empty() {
        let cigar = Cigar::default();
        assert!(cigar.is_empty());

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        assert!(!cigar.is_empty());
    }

    #[test]
    fn test_fmt() {
        let cigar = Cigar::default();
        assert!(cigar.to_string().is_empty());

        let cigar: Cigar = [
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ]
        .into_iter()
        .collect();

        assert_eq!(cigar.to_string(), "1M13N144S");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "1M13N144S".parse::<Cigar>(),
            Ok([
                Op::new(Kind::Match, 1),
                Op::new(Kind::Skip, 13),
                Op::new(Kind::SoftClip, 144),
            ]
            .into_iter()
            .collect())
        );

        assert_eq!("".parse::<Cigar>(), Err(ParseError::Empty));
        assert_eq!("8M13".parse::<Cigar>(), Err(ParseError::Invalid));

        assert!(matches!(
            "*".parse::<Cigar>(),
            Err(ParseError::InvalidOp(_))
        ));
    }
}
