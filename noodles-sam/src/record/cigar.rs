//! SAM CIGAR and operations.

pub mod op;

use std::{io, ops::Deref};

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
    type Error = io::Error;

    fn try_from(ops: Vec<Op>) -> Result<Self, Self::Error> {
        Ok(Self(ops))
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
}
