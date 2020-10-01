//! Tabix index and fields.

pub mod builder;
pub mod header;
mod indexer;
pub mod reference_sequence;

pub use self::{
    builder::Builder, header::Header, indexer::Indexer, reference_sequence::ReferenceSequence,
};

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

    /// Returns a list of reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::{self as tabix, index::ReferenceSequence};
    ///
    /// let reference_sequences = vec![ReferenceSequence::new(Vec::new(), Vec::new())];
    ///
    /// let index = tabix::Index::builder()
    ///     .set_reference_sequences(reference_sequences)
    ///     .build();
    ///
    /// assert_eq!(index.reference_sequences().len(), 1);
    /// ```
    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unmapped reads in the associated file.
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
    pub fn unmapped_read_count(&self) -> Option<u64> {
        self.unmapped_read_count
    }
}

impl Default for Index {
    fn default() -> Self {
        Builder::default().build()
    }
}
