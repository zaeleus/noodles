//! Tabix index header.

mod builder;
pub mod format;

pub use self::{builder::Builder, format::Format};

use indexmap::IndexSet;

/// An ordered set of reference sequence names.
pub type ReferenceSequenceNames = IndexSet<String>;

/// A tabix index header.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    format: Format,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: Option<usize>,
    line_comment_prefix: u8,
    line_skip_count: u32,
    /// A set of ordered reference sequence names.
    pub reference_sequence_names: ReferenceSequenceNames,
}

impl Header {
    /// Creates a tabix index header builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let builder = csi::index::Header::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::{self as csi, index::header::Format};
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_format(Format::Vcf)
    ///     .build();
    ///
    /// assert_eq!(header.format(), Format::Vcf);
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
    /// use noodles_csi as csi;
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_reference_sequence_name_index(1)
    ///     .build();
    ///
    /// assert_eq!(header.reference_sequence_name_index(), 1);
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
    /// use noodles_csi as csi;
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_start_position_index(4)
    ///     .build();
    ///
    /// assert_eq!(header.start_position_index(), 4);
    /// ```
    pub fn start_position_index(&self) -> usize {
        self.start_position_index
    }

    /// Returns the end position field index.
    ///
    /// This index is 1-based. It is not set if the format does not have a column for the end
    /// position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_end_position_index(Some(5))
    ///     .build();
    ///
    /// assert_eq!(header.end_position_index(), Some(5));
    /// ```
    pub fn end_position_index(&self) -> Option<usize> {
        self.end_position_index
    }

    /// Returns the line comment prefix.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_line_comment_prefix(b'#')
    ///     .build();
    ///
    /// assert_eq!(header.line_comment_prefix(), b'#');
    /// ```
    pub fn line_comment_prefix(&self) -> u8 {
        self.line_comment_prefix
    }

    /// Returns the number of lines to skip.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_line_skip_count(0)
    ///     .build();
    ///
    /// assert_eq!(header.line_skip_count(), 0);
    /// ```
    pub fn line_skip_count(&self) -> u32 {
        self.line_skip_count
    }

    /// Returns the reference sequence names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::{self as csi, index::header::ReferenceSequenceNames};
    ///
    /// let reference_sequence_names = ReferenceSequenceNames::new();
    ///
    /// let header = csi::index::Header::builder()
    ///     .set_reference_sequence_names(reference_sequence_names.clone())
    ///     .build();
    ///
    /// assert_eq!(header.reference_sequence_names(), &reference_sequence_names);
    /// ```
    pub fn reference_sequence_names(&self) -> &ReferenceSequenceNames {
        &self.reference_sequence_names
    }
}

impl Default for Header {
    fn default() -> Self {
        Builder::default().build()
    }
}
