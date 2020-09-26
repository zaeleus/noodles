//! Tabix index builder.

use super::{format::CoordinateSystem, Format, Index, ReferenceSequence};

/// A tabix index builder.
pub struct Builder {
    format: Format,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: Option<usize>,
    line_comment_prefix: u8,
    line_skip_count: u32,
    reference_sequence_names: Vec<String>,
    reference_sequences: Vec<ReferenceSequence>,
    unmapped_read_count: Option<u64>,
}

impl Builder {
    /// Creates a builder that targets the GFF format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let builder = tabix::index::builder::Builder::gff();
    /// ```
    pub fn gff() -> Self {
        Builder::default()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_name_index(1)
            .set_start_position_index(4)
            .set_end_position_index(Some(5))
            .set_line_comment_prefix(b'#')
            .set_line_skip_count(0)
    }

    /// Creates a builder that targets the SAM format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let builder = tabix::index::builder::Builder::sam();
    /// ```
    pub fn sam() -> Self {
        Builder::default()
            .set_format(Format::Sam)
            .set_reference_sequence_name_index(3)
            .set_start_position_index(4)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'@')
            .set_line_skip_count(0)
    }

    /// Creates a builder that targets the VCF format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let builder = tabix::index::builder::Builder::vcf();
    /// ```
    pub fn vcf() -> Self {
        Builder::default()
            .set_format(Format::Vcf)
            .set_reference_sequence_name_index(1)
            .set_start_position_index(2)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'#')
            .set_line_skip_count(0)
    }

    /// Sets a tabix index format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::{self as tabix, index::Format};
    /// let index = tabix::Index::builder().set_format(Format::Vcf).build();
    /// assert_eq!(index.format(), Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Sets a tabix index reference sequence name index.
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
    pub fn set_reference_sequence_name_index(
        mut self,
        reference_sequence_name_index: usize,
    ) -> Self {
        self.reference_sequence_name_index = reference_sequence_name_index;
        self
    }

    /// Sets a tabix index start position index.
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
    pub fn set_start_position_index(mut self, start_position_index: usize) -> Self {
        self.start_position_index = start_position_index;
        self
    }

    /// Sets a tabix index end position index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder()
    ///     .set_end_position_index(Some(5))
    ///     .build();
    ///
    /// assert_eq!(index.end_position_index(), Some(5));
    /// ```
    pub fn set_end_position_index(mut self, end_position_index: Option<usize>) -> Self {
        self.end_position_index = end_position_index;
        self
    }

    /// Sets a tabix index line comment prefix.
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
    pub fn set_line_comment_prefix(mut self, line_comment_prefix: u8) -> Self {
        self.line_comment_prefix = line_comment_prefix;
        self
    }

    /// Sets a tabix index line skip count.
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
    /// assert_eq!(index.line_skip_count(), 0);
    /// ```
    pub fn set_line_skip_count(mut self, line_skip_count: u32) -> Self {
        self.line_skip_count = line_skip_count;
        self
    }

    /// Sets tabix index reference sequence names.
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
    pub fn set_reference_sequence_names(mut self, reference_sequence_names: Vec<String>) -> Self {
        self.reference_sequence_names = reference_sequence_names;
        self
    }

    /// Sets tabix index reference sequences.
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
    pub fn set_reference_sequences(mut self, reference_sequences: Vec<ReferenceSequence>) -> Self {
        self.reference_sequences = reference_sequences;
        self
    }

    /// Sets a tabix index unmapped read count.
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
            format: self.format,
            reference_sequence_name_index: self.reference_sequence_name_index,
            start_position_index: self.start_position_index,
            end_position_index: self.end_position_index,
            line_comment_prefix: self.line_comment_prefix,
            line_skip_count: self.line_skip_count,
            reference_sequence_names: self.reference_sequence_names,
            reference_sequences: self.reference_sequences,
            unmapped_read_count: self.unmapped_read_count,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Builder {
            format: Format::Generic(CoordinateSystem::Gff),
            reference_sequence_name_index: 1,
            start_position_index: 4,
            end_position_index: Some(5),
            line_comment_prefix: b'#',
            line_skip_count: 0,
            reference_sequence_names: Vec::new(),
            reference_sequences: Vec::new(),
            unmapped_read_count: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gff() {
        let builder = Builder::gff();
        assert_eq!(builder.format, Format::Generic(CoordinateSystem::Gff));
        assert_eq!(builder.reference_sequence_name_index, 1);
        assert_eq!(builder.start_position_index, 4);
        assert_eq!(builder.end_position_index, Some(5));
        assert_eq!(builder.line_comment_prefix, b'#');
        assert_eq!(builder.line_skip_count, 0);
        assert!(builder.reference_sequence_names.is_empty());
        assert!(builder.reference_sequences.is_empty());
        assert!(builder.unmapped_read_count.is_none());
    }

    #[test]
    fn test_sam() {
        let builder = Builder::sam();
        assert_eq!(builder.format, Format::Sam);
        assert_eq!(builder.reference_sequence_name_index, 3);
        assert_eq!(builder.start_position_index, 4);
        assert_eq!(builder.end_position_index, None);
        assert_eq!(builder.line_comment_prefix, b'@');
        assert_eq!(builder.line_skip_count, 0);
        assert!(builder.reference_sequence_names.is_empty());
        assert!(builder.reference_sequences.is_empty());
        assert!(builder.unmapped_read_count.is_none());
    }

    #[test]
    fn test_vcf() {
        let builder = Builder::vcf();
        assert_eq!(builder.format, Format::Vcf);
        assert_eq!(builder.reference_sequence_name_index, 1);
        assert_eq!(builder.start_position_index, 2);
        assert_eq!(builder.end_position_index, None);
        assert_eq!(builder.line_comment_prefix, b'#');
        assert_eq!(builder.line_skip_count, 0);
        assert!(builder.reference_sequence_names.is_empty());
        assert!(builder.reference_sequences.is_empty());
        assert!(builder.unmapped_read_count.is_none());
    }
}
