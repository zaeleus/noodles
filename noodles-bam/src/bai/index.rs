//! BAM index and fields

pub mod reference;

pub use self::reference::Reference;

/// A BAM index.
#[derive(Debug, Default)]
pub struct Index {
    references: Vec<Reference>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Creates a new BAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let index = bai::Index::new(Vec::new(), None);
    /// ```
    pub fn new(references: Vec<Reference>, n_no_coor: Option<u64>) -> Self {
        Self {
            references,
            n_no_coor,
        }
    }

    /// Returns a list of reference sequences.
    ///
    /// This list is parallel to the reference sequences defined in the associated BAM file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let index = bai::Index::default();
    /// assert!(index.references().is_empty());
    /// ```
    pub fn references(&self) -> &[Reference] {
        &self.references
    }

    /// Returns the number of unmapped reads in the associated BAM file.
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
    pub fn unplaced_unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}
