use noodles_core::Position;
use noodles_gff::feature::{
    record::{Phase, Strand},
    record_buf::Attributes,
};

use super::{RecordBuf, MISSING_FIELD};

/// A GTF record builder.
#[derive(Debug)]
pub struct Builder {
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

impl Builder {
    /// Sets the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().set_reference_sequence_name("sq0").build();
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// ```
    pub fn set_reference_sequence_name<N>(mut self, reference_sequence_name: N) -> Self
    where
        N: Into<String>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the source.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().set_source("NOODLES").build();
    /// assert_eq!(record.source(), "NOODLES");
    /// ```
    pub fn set_source<S>(mut self, source: S) -> Self
    where
        S: Into<String>,
    {
        self.source = source.into();
        self
    }

    /// Sets the feature type.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().set_type("exon").build();
    /// assert_eq!(record.ty(), "exon");
    /// ```
    pub fn set_type<T>(mut self, ty: T) -> Self
    where
        T: Into<String>,
    {
        self.ty = ty.into();
        self
    }

    /// Sets the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    /// let start = Position::MIN;
    /// let record = gtf::RecordBuf::builder().set_start(start).build();
    /// assert_eq!(record.start(), start);
    /// ```
    pub fn set_start(mut self, start: Position) -> Self {
        self.start = start;
        self
    }

    /// Sets the end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    /// let end = Position::MIN;
    /// let record = gtf::RecordBuf::builder().set_end(end).build();
    /// assert_eq!(record.end(), end);
    /// ```
    pub fn set_end(mut self, end: Position) -> Self {
        self.end = end;
        self
    }

    /// Sets the score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().set_score(1.0).build();
    /// assert_eq!(record.score(), Some(1.0));
    /// ```
    pub fn set_score(mut self, score: f32) -> Self {
        self.score = Some(score);
        self
    }

    /// Sets the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::feature::record::Strand;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().set_strand(Strand::Forward).build();
    /// assert_eq!(record.strand(), Strand::Forward);
    /// ```
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Sets the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::feature::record::Phase;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().set_frame(Phase::One).build();
    /// assert_eq!(record.frame(), Some(Phase::One));
    /// ```
    pub fn set_frame(mut self, frame: Phase) -> Self {
        self.frame = Some(frame);
        self
    }

    /// Sets the attributes;
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::feature::record_buf::{attributes::field::Value, Attributes};
    /// use noodles_gtf as gtf;
    ///
    /// let attributes: Attributes = [(String::from("gene_id"), Value::from("g0"))]
    ///     .into_iter()
    ///     .collect();
    ///
    /// let record = gtf::RecordBuf::builder().set_attributes(attributes.clone()).build();
    ///
    /// assert_eq!(record.attributes(), &attributes);
    /// ```
    pub fn set_attributes(mut self, attributes: Attributes) -> Self {
        self.attributes = attributes;
        self
    }

    /// Builds the GTF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::builder().build();
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
            frame: self.frame,
            attributes: self.attributes,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            reference_sequence_name: MISSING_FIELD.into(),
            source: MISSING_FIELD.into(),
            ty: MISSING_FIELD.into(),
            start: Position::MIN,
            end: Position::MIN,
            score: None,
            strand: Strand::None,
            frame: None,
            attributes: Attributes::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert_eq!(builder.reference_sequence_name, ".");
        assert_eq!(builder.source, ".");
        assert_eq!(builder.ty, ".");
        assert_eq!(builder.start, Position::MIN);
        assert_eq!(builder.end, Position::MIN);
        assert!(builder.score.is_none());
        assert_eq!(builder.strand, Strand::None);
        assert!(builder.frame.is_none());
        assert!(builder.attributes.is_empty());
    }
}
