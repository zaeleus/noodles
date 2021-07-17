//! BAM index and fields.

mod builder;
pub mod reference_sequence;

pub use self::{builder::Builder, reference_sequence::ReferenceSequence};

use noodles_csi::BinningIndex;

/// A BAM index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index {
    reference_sequences: Vec<ReferenceSequence>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Creates a BAM index builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let builder = bai::Index::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a BAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let index = bai::Index::new(Vec::new(), None);
    /// ```
    pub fn new(reference_sequences: Vec<ReferenceSequence>, n_no_coor: Option<u64>) -> Self {
        Self {
            reference_sequences,
            n_no_coor,
        }
    }

    /// Returns the number of unplaced unmapped reads in the associated BAM file.
    ///
    /// An unplaced unmapped read is a read that is has neither a reference sequence ID nor
    /// position set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    ///
    /// let index = bai::Index::default();
    /// assert_eq!(index.unplaced_unmapped_read_count(), None);
    ///
    /// let index = bai::Index::new(Vec::new(), Some(13));
    /// assert_eq!(index.unplaced_unmapped_read_count(), Some(13));
    /// ```
    #[deprecated(
        since = "0.2.0",
        note = "Use `unplaced_unmapped_record_count` instead."
    )]
    pub fn unplaced_unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}

impl BinningIndex<ReferenceSequence> for Index {
    /// Returns a list of indexed reference sequences.
    ///
    /// This list is parallel to the reference sequences defined in the associated BAM file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use noodles_csi::BinningIndex;
    /// let index = bai::Index::default();
    /// assert!(index.reference_sequences().is_empty());
    /// ```
    fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unplaced, unmapped records in the associated file.
    ///
    /// An unplaced, unmapped record is one that is has neither a reference sequence ID nor
    /// position set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use noodles_csi::BinningIndex;
    /// let index = bai::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}
