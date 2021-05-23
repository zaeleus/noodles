//! Tabix index reference sequence and fields.

pub mod bin;
mod builder;
pub mod metadata;

pub use self::{bin::Bin, metadata::Metadata};

pub(crate) use self::builder::Builder;

use bit_vec::BitVec;
use noodles_bgzf as bgzf;

const WINDOW_SIZE: i32 = 16384;

/// A tabix index reference sequence.
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

    /// Creates a tabix index reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
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

    /// Returns the list of bins in the reference sequence.
    ///
    /// This list does include the metadata pseudo-bin (bin 37450). Use [`Self::metadata`] instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    /// Returns the list of 16 kbp intervals that make up the linear index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.intervals().is_empty());
    /// ```
    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    /// Returns metadata for this reference sequence.
    ///
    /// Metadata is parsed from the optional pseudo-bin 37450.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::VirtualPosition;
    /// use noodles_tabix::index::{reference_sequence::Metadata, ReferenceSequence};
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
    pub fn metadata(&self) -> Option<&Metadata> {
        self.metadata.as_ref()
    }

    /// Returns a list of bins in this reference sequence that intersects the given range.
    ///
    /// `start` and `end` are 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// let query_bins = reference_sequence.query(8, 13);
    /// assert!(query_bins.is_empty());
    /// ```
    pub fn query(&self, start: i32, end: i32) -> Vec<&Bin> {
        let region_bins = region_to_bins((start - 1) as usize, end as usize);

        self.bins()
            .iter()
            .filter(|b| region_bins[b.id() as usize])
            .collect()
    }

    /// Finds in minimum start offset in the linear index for a given start position.
    ///
    /// `start` is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert_eq!(reference_sequence.min_offset(13), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn min_offset(&self, start: i32) -> bgzf::VirtualPosition {
        let i = ((start - 1) / WINDOW_SIZE) as usize;
        self.intervals.get(i).copied().unwrap_or_default()
    }
}

// 0-based, [start, end)
fn region_to_bins(start: usize, mut end: usize) -> BitVec {
    end -= 1;

    let mut bins = BitVec::from_elem(bin::MAX_ID, false);
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
        let mut expected = BitVec::from_elem(bin::MAX_ID, false);
        for &k in &[0, 1, 9, 73, 585, 4681] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);

        // [63245986, 63245986]
        let actual = region_to_bins(63245985, 63255986);
        let mut expected = BitVec::from_elem(bin::MAX_ID, false);
        for &k in &[0, 1, 16, 133, 1067, 8541] {
            expected.set(k, true);
        }
        assert_eq!(actual, expected);
    }
}
