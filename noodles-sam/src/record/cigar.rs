//! SAM CIGAR and operations.

use std::io;

use crate::alignment::record::cigar::Op;

/// A SAM record CIGAR.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Cigar(Vec<Op>);

impl Cigar {
    /// Calculates the alignment span over the reference sequence.
    ///
    /// This sums the lengths of the CIGAR operations that consume the reference sequence, i.e.,
    /// alignment matches (`M`), deletions from the reference (`D`), skipped reference regions
    /// (`S`), sequence matches (`=`), and sequence mismatches (`X`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     alignment::record::cigar::{op::Kind, Op},
    ///     record::Cigar,
    /// };
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
    /// use noodles_sam::{
    ///     alignment::record::cigar::{op::Kind, Op},
    ///     record::Cigar,
    /// };

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
}

impl crate::alignment::record::Cigar for Cigar {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Op>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}

impl crate::alignment::record::Cigar for &Cigar {
    fn is_empty(&self) -> bool {
        (*self).is_empty()
    }

    fn len(&self) -> usize {
        (*self).len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Op>> + '_> {
        (*self).iter()
    }
}

impl AsRef<[Op]> for Cigar {
    fn as_ref(&self) -> &[Op] {
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

impl From<Vec<Op>> for Cigar {
    fn from(ops: Vec<Op>) -> Self {
        Self(ops)
    }
}

impl From<Cigar> for Vec<Op> {
    fn from(cigar: Cigar) -> Self {
        cigar.0
    }
}
