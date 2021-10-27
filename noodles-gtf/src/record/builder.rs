use super::{Attributes, Frame, Record, Strand, NULL_FIELD};

/// A GTF record builder.
#[derive(Debug)]
pub struct Builder {
    reference_sequence_name: String,
    source: String,
    ty: String,
    start: i32,
    end: i32,
    score: Option<f32>,
    strand: Option<Strand>,
    frame: Option<Frame>,
    attributes: Attributes,
}

impl Builder {
    /// Sets the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::Record::builder().set_reference_sequence_name("sq0").build();
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
    /// let record = gtf::Record::builder().set_source("NOODLES").build();
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
    /// let record = gtf::Record::builder().set_type("exon").build();
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
    /// use noodles_gtf as gtf;
    /// let record = gtf::Record::builder().set_start(8).build();
    /// assert_eq!(record.start(), 8);
    /// ```
    pub fn set_start(mut self, start: i32) -> Self {
        self.start = start;
        self
    }

    /// Sets the end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::Record::builder().set_end(13).build();
    /// assert_eq!(record.end(), 13);
    /// ```
    pub fn set_end(mut self, end: i32) -> Self {
        self.end = end;
        self
    }

    /// Sets the score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::Record::builder().set_score(1.0).build();
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
    /// use noodles_gtf::{self as gtf, record::Strand};
    /// let record = gtf::Record::builder().set_strand(Strand::Forward).build();
    /// assert_eq!(record.strand(), Some(Strand::Forward));
    /// ```
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    /// Sets the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::{self as gtf, record::Frame};
    /// let frame = Frame::try_from(0)?;
    /// let record = gtf::Record::builder().set_frame(frame).build();
    /// assert_eq!(record.frame(), Some(frame));
    /// Ok::<_, gtf::record::frame::ParseError>(())
    /// ```
    pub fn set_frame(mut self, frame: Frame) -> Self {
        self.frame = Some(frame);
        self
    }

    /// Sets the attributes;
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::{self as gtf, record::{attributes::Entry, Attributes}};
    /// let attributes = Attributes::from(vec![Entry::new("gene_id", "g0")]);
    /// let record = gtf::Record::builder().set_attributes(attributes.clone()).build();
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
    /// let record = gtf::Record::builder().build();
    /// ```
    pub fn build(self) -> Record {
        Record {
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
            reference_sequence_name: NULL_FIELD.into(),
            source: NULL_FIELD.into(),
            ty: NULL_FIELD.into(),
            start: 1,
            end: 1,
            score: None,
            strand: None,
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
        assert_eq!(builder.start, 1);
        assert_eq!(builder.end, 1);
        assert!(builder.score.is_none());
        assert!(builder.strand.is_none());
        assert!(builder.frame.is_none());
        assert!(builder.attributes.is_empty());
    }
}
