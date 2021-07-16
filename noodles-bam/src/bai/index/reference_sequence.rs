//! BAM index reference sequence and fields.

pub mod bin;
mod builder;

pub(crate) use self::builder::Builder;

pub use self::bin::Bin;

use std::{
    error, fmt,
    ops::{Bound, RangeBounds},
};

use bit_vec::BitVec;
use noodles_bgzf as bgzf;
use noodles_csi::{index::reference_sequence::Metadata, BinningIndexReferenceSequence};

const MIN_SHIFT: i32 = 14;
const DEPTH: i32 = 5;
const MIN_POSITION: i32 = 1;
const MAX_POSITION: i32 = 1 << (MIN_SHIFT + 3 * DEPTH);

const WINDOW_SIZE: i32 = 1 << MIN_SHIFT;

/// An error returned when a query fails.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum QueryError {
    /// The start position is invalid.
    InvalidStartPosition(i32),
    /// The end position is invalid.
    InvalidEndPosition(i32),
}

impl error::Error for QueryError {}

impl fmt::Display for QueryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStartPosition(start) => {
                write!(
                    f,
                    "expected start position >= {}, got {}",
                    MIN_POSITION, start
                )
            }
            Self::InvalidEndPosition(end) => {
                write!(f, "expected end position <= {}, got {}", MAX_POSITION, end)
            }
        }
    }
}

/// A reference sequence in the BAM index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
    metadata: Option<Metadata>,
}

impl ReferenceSequence {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a BAM index reference seqeuence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// ```
    pub fn new(
        bins: Vec<Bin>,
        intervals: Vec<bgzf::VirtualPosition>,
        metadata: Option<Metadata>,
    ) -> Self {
        Self {
            bins,
            intervals,
            metadata,
        }
    }

    /// Returns the list of bins in this reference sequence.
    ///
    /// This list does not include the metadata pseudo-bin (bin 37450). Use [`Self::metadata`]
    /// instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    /// Returns a list of 16 kbp intervals that make up the linear index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.intervals().is_empty());
    /// ```
    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    /// Returns a list of bins in this reference seqeunce that intersect the given range.
    ///
    /// The interval values are 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_bam::bai::index::reference_sequence;
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// let query_bins = reference_sequence.query(8..=13)?;
    /// assert!(query_bins.is_empty());
    /// # Ok::<(), reference_sequence::QueryError>(())
    /// ```
    pub fn query<B>(&self, interval: B) -> Result<Vec<&Bin>, QueryError>
    where
        B: RangeBounds<i32>,
    {
        let start = match interval.start_bound() {
            Bound::Included(s) => *s,
            Bound::Excluded(s) => *s + 1,
            Bound::Unbounded => MIN_POSITION,
        };

        if start < MIN_POSITION {
            return Err(QueryError::InvalidStartPosition(start));
        }

        let end = match interval.end_bound() {
            Bound::Included(e) => *e,
            Bound::Excluded(e) => *e - 1,
            Bound::Unbounded => MAX_POSITION,
        };

        if end > MAX_POSITION {
            return Err(QueryError::InvalidEndPosition(end));
        }

        let region_bins = region_to_bins((start - 1) as usize, end as usize);

        let query_bins = self
            .bins()
            .iter()
            .filter(|b| region_bins[b.id() as usize])
            .collect();

        Ok(query_bins)
    }

    /// Finds in minimum start offset in the linear index for a given start position.
    ///
    /// `start` is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert_eq!(reference_sequence.min_offset(13), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn min_offset(&self, start: i32) -> bgzf::VirtualPosition {
        let i = ((start - 1) / WINDOW_SIZE) as usize;
        self.intervals.get(i).copied().unwrap_or_default()
    }
}

impl BinningIndexReferenceSequence for ReferenceSequence {
    /// Returns the optional metadata for the reference sequence.
    ///
    /// Metadata is parsed from the optional pseudo-bin 37450.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::VirtualPosition;
    /// use noodles_csi::{index::reference_sequence::Metadata, BinningIndexReferenceSequence};
    /// use noodles_bam::bai::index::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.metadata().is_none());
    ///
    /// let reference_sequence = ReferenceSequence::new(
    ///     Vec::new(),
    ///     Vec::new(),
    ///     Some(Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0))
    /// );
    /// assert!(reference_sequence.metadata().is_some());
    /// ```
    fn metadata(&self) -> Option<&Metadata> {
        self.metadata.as_ref()
    }
}

// 0-based, [start, end)
fn region_to_bins(start: usize, mut end: usize) -> BitVec {
    end -= 1;

    let mut bins = BitVec::from_elem(bin::MAX_ID as usize, false);
    bins.set(0, true);

    for k in (1 + (start >> 26))..=(1 + (end >> 26)) {
        bins.set(k, true);
    }

    for k in (9 + (start >> 23))..=(9 + (end >> 23)) {
        bins.set(k, true);
    }

    for k in (73 + (start >> 20))..=(73 + (end >> 20)) {
        bins.set(k, true);
    }

    for k in (585 + (start >> 17))..=(585 + (end >> 17)) {
        bins.set(k, true);
    }

    for k in (4681 + (start >> 14))..=(4681 + (end >> 14)) {
        bins.set(k, true);
    }

    bins
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query() {
        let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);

        assert_eq!(
            reference_sequence.query(0..=8),
            Err(QueryError::InvalidStartPosition(0))
        );

        assert_eq!(
            reference_sequence.query(1..=i32::MAX),
            Err(QueryError::InvalidEndPosition(i32::MAX))
        );
    }

    #[test]
    fn test_region_to_bins() {
        // [8, 13]
        let actual = region_to_bins(7, 13);
        let mut expected = BitVec::from_elem(bin::MAX_ID as usize, false);
        for &k in &[0, 1, 9, 73, 585, 4681] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);

        // [63245986, 63245986]
        let actual = region_to_bins(63245985, 63255986);
        let mut expected = BitVec::from_elem(bin::MAX_ID as usize, false);
        for &k in &[0, 1, 16, 133, 1067, 8541] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);
    }
}
