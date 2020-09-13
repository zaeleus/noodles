//! BAM index reference sequence and fields.

pub mod bin;
mod builder;
pub mod metadata;

pub(crate) use self::builder::Builder;

pub use self::{bin::Bin, metadata::Metadata};

use std::convert::TryFrom;

use bit_vec::BitVec;
use noodles_bgzf as bgzf;

use self::bin::MAX_BIN;

const WINDOW_SIZE: u64 = 16384;

/// A reference sequence in the BAM index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
}

impl ReferenceSequence {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a new BAM index reference seqeuence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// ```
    pub fn new(bins: Vec<Bin>, intervals: Vec<bgzf::VirtualPosition>) -> Self {
        Self { bins, intervals }
    }

    /// Returns the list of bins in this reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
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
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// assert!(reference_sequence.intervals().is_empty());
    /// ```
    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    /// Returns metadata for this reference sequence.
    ///
    /// Metadata is parsed from optional bin 37450.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::{
    ///     reference_sequence::{bin::Chunk, Bin, Metadata},
    ///     ReferenceSequence,
    /// };
    /// use noodles_bgzf::VirtualPosition;
    ///
    /// let reference_sequences = ReferenceSequence::new(
    ///     vec![
    ///         Bin::new(0, Vec::new()),
    ///         Bin::new(37450, vec![
    ///             Chunk::new(VirtualPosition::from(610), VirtualPosition::from(1597)),
    ///             Chunk::new(VirtualPosition::from(55), VirtualPosition::from(0)),
    ///         ]),
    ///     ],
    ///     Vec::new(),
    /// );
    ///
    /// assert_eq!(reference_sequences.metadata(), Some(Metadata::new(
    ///     VirtualPosition::from(610),
    ///     VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// )));
    /// ```
    pub fn metadata(&self) -> Option<Metadata> {
        self.bins()
            .iter()
            .find(|b| b.id() == metadata::MAGIC_NUMBER)
            .and_then(|bin| Metadata::try_from(bin).ok())
    }

    /// Returns a list of bins in this reference seqeunce that intersect the given range.
    ///
    /// `start` and `end` are 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// let query_bins = reference_sequence.query(8, 13);
    /// assert!(query_bins.is_empty());
    /// ```
    pub fn query(&self, start: u64, end: u64) -> Vec<&Bin> {
        let region_bins = region_to_bins((start - 1) as usize, end as usize);

        let mut query_bins = Vec::new();

        for bin in self.bins() {
            let bin_id = bin.id() as usize;

            // Only accept bin numbers [0, MAX_BIN), which skips the psuedo-bin at 37450.
            if bin_id < region_bins.len() && region_bins[bin_id] {
                query_bins.push(bin);
            }
        }

        query_bins
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
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// assert_eq!(reference_sequence.min_offset(13), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn min_offset(&self, start: u64) -> bgzf::VirtualPosition {
        let i = ((start - 1) / WINDOW_SIZE) as usize;
        self.intervals.get(i).copied().unwrap_or_default()
    }
}

// 0-based, [start, end)
fn region_to_bins(start: usize, mut end: usize) -> BitVec {
    end -= 1;

    let mut bins = BitVec::from_elem(MAX_BIN, false);
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
    fn test_region_to_bins() {
        // [8, 13]
        let actual = region_to_bins(7, 13);
        let mut expected = BitVec::from_elem(MAX_BIN, false);
        for &k in &[0, 1, 9, 73, 585, 4681] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);

        // [63245986, 63245986]
        let actual = region_to_bins(63245985, 63255986);
        let mut expected = BitVec::from_elem(MAX_BIN, false);
        for &k in &[0, 1, 16, 133, 1067, 8541] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);
    }
}
