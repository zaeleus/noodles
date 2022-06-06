//! BAM index reference sequence and fields.

pub mod bin;
mod builder;

pub(crate) use self::builder::Builder;

pub use self::bin::Bin;

use std::io;

use bit_vec::BitVec;
use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Position};
use noodles_csi::{binning_index::ReferenceSequenceExt, index::reference_sequence::Metadata};

use super::{resolve_interval, MIN_SHIFT};

const WINDOW_SIZE: usize = 1 << MIN_SHIFT;

/// A reference sequence in the BAM index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
    metadata: Option<Metadata>,
}

impl ReferenceSequence {
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

    /// Returns a list of bins in this reference sequence that intersect the given range.
    ///
    /// The interval values are 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// use noodles_core::Position;
    ///
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    ///
    /// let query_bins = reference_sequence.query(start..=end)?;
    /// assert!(query_bins.is_empty());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<I>(&self, interval: I) -> io::Result<Vec<&Bin>>
    where
        I: Into<Interval>,
    {
        let (start, end) = resolve_interval(interval)?;
        let region_bins = region_to_bins(start, end);

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
    /// use noodles_core::Position;
    ///
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// let start = Position::try_from(13)?;
    ///
    /// assert_eq!(
    ///     reference_sequence.min_offset(start),
    ///     bgzf::VirtualPosition::from(0)
    /// );
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn min_offset(&self, start: Position) -> bgzf::VirtualPosition {
        let i = (usize::from(start) - 1) / WINDOW_SIZE;
        self.intervals.get(i).copied().unwrap_or_default()
    }
}

impl ReferenceSequenceExt for ReferenceSequence {
    /// Returns the optional metadata for the reference sequence.
    ///
    /// Metadata is parsed from the optional pseudo-bin 37450.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::VirtualPosition;
    /// use noodles_csi::{binning_index::ReferenceSequenceExt, index::reference_sequence::Metadata};
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

    /// Returns the start position of the first record in the last linear bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::binning_index::ReferenceSequenceExt;
    /// use noodles_bam::bai::index::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::default();
    /// assert!(reference_sequence.first_record_in_last_linear_bin_start_position().is_none());
    ///
    /// let intervals = vec![bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13)];
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), intervals, None);
    /// assert_eq!(
    ///     reference_sequence.first_record_in_last_linear_bin_start_position(),
    ///     Some(bgzf::VirtualPosition::from(13))
    /// );
    /// ```
    fn first_record_in_last_linear_bin_start_position(&self) -> Option<bgzf::VirtualPosition> {
        self.intervals().last().copied()
    }
}

fn region_to_bins(start: Position, end: Position) -> BitVec {
    // 0-based, [start, end)
    let start = usize::from(start) - 1;
    let end = usize::from(end) - 1;

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

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_query() -> Result<(), noodles_core::position::TryFromIntError> {
        let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
        let end = usize::try_from(i32::MAX).and_then(Position::try_from)?;

        assert!(matches!(
            reference_sequence.query(..=end),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_region_to_bins() -> Result<(), noodles_core::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        let actual = region_to_bins(start, end);
        let mut expected = BitVec::from_elem(bin::MAX_ID as usize, false);
        for &k in &[0, 1, 9, 73, 585, 4681] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);

        let start = Position::try_from(63245985)?;
        let end = Position::try_from(63255986)?;
        let actual = region_to_bins(start, end);
        let mut expected = BitVec::from_elem(bin::MAX_ID as usize, false);
        for &k in &[0, 1, 16, 133, 1067, 8541] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);

        Ok(())
    }
}
