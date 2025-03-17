//! GTF record and fields.

mod attributes;
mod builder;
mod convert;

use noodles_core::Position;
use noodles_gff::feature::record::{Phase, Strand};

pub use self::{attributes::Attributes, builder::Builder};

pub(crate) const MISSING_FIELD: &str = ".";

/// A GTF record buffer.
#[derive(Clone, Debug, PartialEq)]
pub struct RecordBuf {
    reference_sequence_name: String,
    source: String,
    ty: String,
    start: Position,
    end: Position,
    score: Option<f32>,
    strand: Strand,
    frame: Option<Phase>,
    attributes: Attributes,
}

impl RecordBuf {
    /// Returns a record builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let builder = gtf::RecordBuf::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the reference sequence name.
    ///
    /// This is also called the "seqname".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.reference_sequence_name(), ".");
    /// ```
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    /// Returns the source.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.source(), ".");
    /// ```
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Returns the feature type.
    ///
    /// This is also simply called "feature".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.ty(), ".");
    /// ```
    pub fn ty(&self) -> &str {
        &self.ty
    }

    /// Returns the start position.
    ///
    /// This position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.start(), Position::MIN);
    /// ```
    pub fn start(&self) -> Position {
        self.start
    }

    /// Returns the end position.
    ///
    /// This position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.end(), Position::MIN);
    /// ```
    pub fn end(&self) -> Position {
        self.end
    }

    /// Returns the confidence score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.score().is_none());
    /// ```
    pub fn score(&self) -> Option<f32> {
        self.score
    }

    /// Returns the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::feature::record::Strand;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.strand(), Strand::None);
    /// ```
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.frame().is_none());
    /// ```
    pub fn frame(&self) -> Option<Phase> {
        self.frame
    }

    /// Returns the attributes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.attributes().is_empty());
    /// ```
    pub fn attributes(&self) -> &Attributes {
        &self.attributes
    }
}

impl Default for RecordBuf {
    fn default() -> Self {
        Self::builder().build()
    }
}
