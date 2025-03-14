use noodles_core::Position;

use super::{Attributes, RecordBuf};
use crate::record::{Phase, Strand};

/// A GFF record builder.
#[derive(Debug)]
pub struct Builder {
    reference_sequence_name: String,
    source: String,
    ty: String,
    start: Position,
    end: Position,
    score: Option<f32>,
    strand: Strand,
    phase: Option<Phase>,
    attributes: Attributes,
}

impl Builder {
    /// Creates a GFF record builder.
    ///
    /// Typically, [`RecordBuf::builder`] is used instead of calling this.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let builder = gff::feature::RecordBuf::builder();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets a GFF record reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    ///
    /// let record = gff::feature::RecordBuf::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// ```
    pub fn set_reference_sequence_name<N>(mut self, reference_sequence_name: N) -> Self
    where
        N: Into<String>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets a GFF record source.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    ///
    /// let record = gff::feature::RecordBuf::builder()
    ///     .set_source("NOODLES")
    ///     .build();
    ///
    /// assert_eq!(record.source(), "NOODLES");
    /// ```
    pub fn set_source<S>(mut self, source: S) -> Self
    where
        S: Into<String>,
    {
        self.source = source.into();
        self
    }

    /// Sets a GFF record feature type.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    ///
    /// let record = gff::feature::RecordBuf::builder()
    ///     .set_type("gene")
    ///     .build();
    ///
    /// assert_eq!(record.ty(), "gene");
    /// ```
    pub fn set_type<T>(mut self, ty: T) -> Self
    where
        T: Into<String>,
    {
        self.ty = ty.into();
        self
    }

    /// Sets a GFF record start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gff as gff;
    /// let start = Position::MIN;
    /// let record = gff::feature::RecordBuf::builder().set_start(start).build();
    /// assert_eq!(record.start(), start);
    /// ```
    pub fn set_start(mut self, start: Position) -> Self {
        self.start = start;
        self
    }

    /// Sets a GFF record end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gff as gff;
    /// let end = Position::MIN;
    /// let record = gff::feature::RecordBuf::builder().set_end(end).build();
    /// assert_eq!(record.end(), end);
    /// ```
    pub fn set_end(mut self, end: Position) -> Self {
        self.end = end;
        self
    }

    /// Sets a GFF record score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::feature::RecordBuf::builder().set_score(21.0).build();
    /// assert_eq!(record.score(), Some(21.0));
    /// ```
    pub fn set_score(mut self, score: f32) -> Self {
        self.score = Some(score);
        self
    }

    /// Sets a GFF record strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::{self as gff, record::Strand};
    ///
    /// let record = gff::feature::RecordBuf::builder()
    ///     .set_strand(Strand::Forward)
    ///     .build();
    ///
    /// assert_eq!(record.strand(), Strand::Forward);
    /// ```
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Sets a GFF record phase.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::{self as gff, record::Phase};
    /// let record = gff::feature::RecordBuf::builder().set_phase(Phase::Zero).build();
    /// assert_eq!(record.phase(), Some(Phase::Zero));
    /// ```
    pub fn set_phase(mut self, phase: Phase) -> Self {
        self.phase = Some(phase);
        self
    }

    /// Sets GFF record attributes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::{
    ///     self as gff,
    ///     feature::record_buf::{
    ///         attributes::field::{Tag, Value},
    ///         Attributes,
    ///     },
    /// };
    ///
    /// let attributes: Attributes = [(Tag::from("gene_id"), Value::from("ndls0"))]
    ///     .into_iter()
    ///     .collect();
    ///
    /// let record = gff::feature::RecordBuf::builder()
    ///     .set_attributes(attributes.clone())
    ///     .build();
    ///
    /// assert_eq!(record.attributes(), &attributes);
    /// ```
    pub fn set_attributes(mut self, attributes: Attributes) -> Self {
        self.attributes = attributes;
        self
    }

    /// Builds a GFF record.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let record = gff::feature::RecordBuf::builder().build();
    /// ```
    pub fn build(self) -> RecordBuf {
        RecordBuf {
            reference_sequence_name: self.reference_sequence_name,
            source: self.source,
            ty: self.ty,
            start: self.start,
            end: self.end,
            score: self.score,
            strand: self.strand,
            phase: self.phase,
            attributes: self.attributes,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        const MISSING: &str = ".";

        Self {
            reference_sequence_name: MISSING.into(),
            source: MISSING.into(),
            ty: MISSING.into(),
            start: Position::MIN,
            end: Position::MIN,
            score: None,
            strand: Strand::default(),
            phase: None,
            attributes: Attributes::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let record = Builder::default();

        assert_eq!(record.reference_sequence_name, ".");
        assert_eq!(record.source, ".");
        assert_eq!(record.ty, ".");
        assert_eq!(record.start, Position::MIN);
        assert_eq!(record.end, Position::MIN);
        assert!(record.score.is_none());
        assert_eq!(record.strand, Strand::default());
        assert!(record.phase.is_none());
        assert!(record.attributes.is_empty());
    }

    #[test]
    fn test_build() -> Result<(), noodles_core::position::TryFromIntError> {
        use crate::feature::record_buf::attributes::field::{Tag, Value};

        let attributes: Attributes = [(Tag::from("gene_id"), Value::from("ndls0"))]
            .into_iter()
            .collect();

        let record = Builder::new()
            .set_reference_sequence_name("sq0")
            .set_source("NOODLES")
            .set_type("CDS")
            .set_start(Position::try_from(8)?)
            .set_end(Position::try_from(13)?)
            .set_score(21.0)
            .set_strand(Strand::Forward)
            .set_phase(Phase::Zero)
            .set_attributes(attributes.clone())
            .build();

        assert_eq!(record.reference_sequence_name(), "sq0");
        assert_eq!(record.source(), "NOODLES");
        assert_eq!(record.ty(), "CDS");
        assert_eq!(record.start(), Position::try_from(8)?);
        assert_eq!(record.end(), Position::try_from(13)?);
        assert_eq!(record.score(), Some(21.0));
        assert_eq!(record.strand(), Strand::Forward);
        assert_eq!(record.phase(), Some(Phase::Zero));
        assert_eq!(record.attributes(), &attributes);

        Ok(())
    }
}
