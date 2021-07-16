//! Coordinate-sorted index and fields.

mod builder;
pub mod reference_sequence;

pub use self::{builder::Builder, reference_sequence::ReferenceSequence};

use super::BinningIndex;

/// A coordinate-sorted index (CSI).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Index {
    min_shift: i32,
    depth: i32,
    aux: Vec<u8>,
    reference_sequences: Vec<ReferenceSequence>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Returns a builder to create an index from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let builder = csi::Index::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the number of bits for the minimum interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert_eq!(index.min_shift(), 14);
    /// ```
    pub fn min_shift(&self) -> i32 {
        self.min_shift
    }

    /// Returns the depth of the binning index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert_eq!(index.depth(), 5);
    /// ```
    pub fn depth(&self) -> i32 {
        self.depth
    }

    /// Returns the auxiliary data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.aux().is_empty());
    /// ```
    pub fn aux(&self) -> &[u8] {
        &self.aux
    }

    /// Returns the number of unmapped records in the associated file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.unmapped_read_count().is_none());
    /// ```
    #[deprecated(
        since = "0.2.0",
        note = "Use `unplaced_unmapped_record_count` instead."
    )]
    pub fn unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}

impl BinningIndex<ReferenceSequence> for Index {
    /// Returns a list of indexed reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::{self as csi, BinningIndex};
    /// let index = csi::Index::default();
    /// assert!(index.reference_sequences().is_empty());
    /// ```
    fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unplaced, unmapped records in the associated file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::{self as csi, BinningIndex};
    /// let index = csi::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}

impl Default for Index {
    fn default() -> Self {
        Self::builder().build()
    }
}
