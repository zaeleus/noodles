//! Tabix index and fields.

mod builder;
pub mod format;
mod indexer;
pub mod reference_sequence;

pub use self::{
    builder::Builder, format::Format, indexer::Indexer, reference_sequence::ReferenceSequence,
};

/// A tabix index.
#[derive(Debug)]
pub struct Index {
    format: Format,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: usize,
    line_comment_prefix: u8,
    line_skip_count: u32,
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

    /// Returns the format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::{self as tabix, index::Format};
    /// let index = tabix::Index::builder().set_format(Format::Vcf).build();
    /// assert_eq!(index.format(), Format::Vcf);
    /// ```
    pub fn format(&self) -> Format {
        self.format
    }

    /// Returns the reference sequence name field index.
    ///
    /// This index is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_reference_sequence_name_index(1)
    ///     .build();
    ///
    /// assert_eq!(index.reference_sequence_name_index(), 1);
    /// ```
    pub fn reference_sequence_name_index(&self) -> usize {
        self.reference_sequence_name_index
    }

    /// Returns the start position field index.
    ///
    /// This index is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_start_position_index(4)
    ///     .build();
    ///
    /// assert_eq!(index.start_position_index(), 4);
    /// ```
    pub fn start_position_index(&self) -> usize {
        self.start_position_index
    }

    /// Returns the end position field index.
    ///
    /// This index is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_end_position_index(5)
    ///     .build();
    ///
    /// assert_eq!(index.end_position_index(), 5);
    /// ```
    pub fn end_position_index(&self) -> usize {
        self.end_position_index
    }

    /// Returns the line comment prefix.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_line_comment_prefix(b'#')
    ///     .build();
    ///
    /// assert_eq!(index.line_comment_prefix(), b'#');
    /// ```
    pub fn line_comment_prefix(&self) -> u8 {
        self.line_comment_prefix
    }

    /// Returns the number of lines to skip.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_line_skip_count(0)
    ///     .build();
    ///
    /// assert_eq!(index.line_skip_count(),0);
    /// ```
    pub fn line_skip_count(&self) -> u32 {
        self.line_skip_count
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
