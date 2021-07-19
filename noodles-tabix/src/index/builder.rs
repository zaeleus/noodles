//! Tabix index builder.

use indexmap::IndexSet;

use super::{Header, Index, ReferenceSequence, ReferenceSequenceNames};

/// A tabix index builder.
pub struct Builder {
    header: Header,
    reference_sequence_names: ReferenceSequenceNames,
    reference_sequences: Vec<ReferenceSequence>,
    unmapped_read_count: Option<u64>,
}

impl Builder {
    /// Sets a header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let header = tabix::index::Header::default();
    /// let index = tabix::Index::builder().set_header(header.clone()).build();
    /// assert_eq!(index.header(), &header);
    /// ```
    pub fn set_header(mut self, header: Header) -> Self {
        self.header = header;
        self
    }

    /// Sets reference sequence names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::{self as tabix, index::ReferenceSequenceNames};
    ///
    /// let reference_sequence_names: ReferenceSequenceNames = vec![String::from("sq0")]
    ///     .into_iter()
    ///     .collect();
    ///
    /// let index = tabix::Index::builder()
    ///     .set_reference_sequence_names(reference_sequence_names.clone())
    ///     .build();
    ///
    /// assert_eq!(index.reference_sequence_names(), &reference_sequence_names);
    /// ```
    pub fn set_reference_sequence_names(
        mut self,
        reference_sequence_names: ReferenceSequenceNames,
    ) -> Self {
        self.reference_sequence_names = reference_sequence_names;
        self
    }

    /// Sets reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::BinningIndex;
    /// use noodles_tabix::{self as tabix, index::ReferenceSequence};
    ///
    /// let reference_sequences = vec![ReferenceSequence::new(Vec::new(), Vec::new(), None)];
    ///
    /// let index = tabix::Index::builder()
    ///     .set_reference_sequences(reference_sequences)
    ///     .build();
    ///
    /// assert_eq!(index.reference_sequences().len(), 1);
    /// ```
    pub fn set_reference_sequences(mut self, reference_sequences: Vec<ReferenceSequence>) -> Self {
        self.reference_sequences = reference_sequences;
        self
    }

    /// Sets an unmapped read count.
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
    pub fn set_unmapped_read_count(mut self, unmapped_read_count: u64) -> Self {
        self.unmapped_read_count = Some(unmapped_read_count);
        self
    }

    /// Builds a tabix index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let index = tabix::Index::builder().build();
    /// ```
    pub fn build(self) -> Index {
        Index {
            header: self.header,
            reference_sequence_names: self.reference_sequence_names,
            reference_sequences: self.reference_sequences,
            unmapped_read_count: self.unmapped_read_count,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            header: Header::builder().build(),
            reference_sequence_names: IndexSet::new(),
            reference_sequences: Vec::new(),
            unmapped_read_count: None,
        }
    }
}
