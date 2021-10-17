use super::{Attributes, Record, Strand};

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
    frame: Option<String>,
    attributes: Attributes,
}

impl Builder {
    /// Sets the reference sequence name.
    pub fn set_reference_sequence_name<N>(mut self, reference_sequence_name: N) -> Self
    where
        N: Into<String>,
    {
        self.reference_sequence_name = reference_sequence_name.into();
        self
    }

    /// Sets the source.
    pub fn set_source<S>(mut self, source: S) -> Self
    where
        S: Into<String>,
    {
        self.source = source.into();
        self
    }

    /// Sets the feature type.
    pub fn set_type<T>(mut self, ty: T) -> Self
    where
        T: Into<String>,
    {
        self.ty = ty.into();
        self
    }

    /// Sets the start position.
    pub fn set_start(mut self, start: i32) -> Self {
        self.start = start;
        self
    }

    /// Sets the end position.
    pub fn set_end(mut self, end: i32) -> Self {
        self.end = end;
        self
    }

    /// Sets the score.
    pub fn set_score(mut self, score: f32) -> Self {
        self.score = Some(score);
        self
    }

    /// Sets the strand.
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    /// Sets the frame.
    pub fn set_frame<F>(mut self, frame: F) -> Self
    where
        F: Into<String>,
    {
        self.frame = Some(frame.into());
        self
    }

    /// Sets the attributes;
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
            reference_sequence_name: String::from("."),
            source: String::from("."),
            ty: String::from("."),
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
