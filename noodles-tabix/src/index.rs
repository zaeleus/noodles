//! Tabix index and fields.

pub mod builder;
pub mod header;
mod indexer;
pub mod reference_sequence;

pub use self::{
    builder::Builder, header::Header, indexer::Indexer, reference_sequence::ReferenceSequence,
};

use noodles_csi::BinningIndex;

/// A tabix index.
#[derive(Debug)]
pub struct Index {
    header: Header,
    reference_sequence_names: Vec<String>,
    reference_sequences: Vec<ReferenceSequence>,
    unmapped_read_count: Option<u64>,
}

impl Index {
    /// Returns a builder to create an index from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let builder = tabix::Index::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns an indexer to create an index from records.
    pub fn indexer() -> Indexer {
        Indexer::default()
    }

    /// Returns the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let header = tabix::index::Header::default();
    /// let index = tabix::Index::builder().set_header(header.clone()).build();
    /// assert_eq!(index.header(), &header);
    /// ```
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Returns the reference sequence names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_reference_sequence_names(vec![String::from("sq0")])
    ///     .build();
    ///
    /// assert_eq!(index.reference_sequence_names(), [String::from("sq0")]);
    /// ```
    pub fn reference_sequence_names(&self) -> &[String] {
        &self.reference_sequence_names
    }

    /// Returns the number of unmapped records in the associated file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_unmapped_read_count(21)
    ///     .build();
    ///
    /// assert_eq!(index.unmapped_read_count(), Some(21));
    /// ```
    #[deprecated(
        since = "0.2.0",
        note = "Use `unplaced_unmapped_record_count` instead."
    )]
    pub fn unmapped_read_count(&self) -> Option<u64> {
        self.unmapped_read_count
    }
}

impl BinningIndex<ReferenceSequence> for Index {
    /// Returns a list of indexed reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::BinningIndex;
    /// use noodles_tabix as tabix;
    /// let index = tabix::Index::default();
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
    /// use noodles_csi::BinningIndex;
    /// use noodles_tabix as tabix;
    /// let index = tabix::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.unmapped_read_count
    }
}

impl Default for Index {
    fn default() -> Self {
        Builder::default().build()
    }
}
