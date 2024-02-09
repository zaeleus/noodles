//! VCF record builder.

use super::{AlternateBases, Filters, Ids, Info, Position, RecordBuf, Samples};

/// A VCF record builder.
#[derive(Debug, PartialEq)]
pub struct Builder {
    chromosome: String,
    position: Position,
    ids: Ids,
    reference_bases: String,
    alternate_bases: AlternateBases,
    quality_score: Option<f32>,
    filters: Option<Filters>,
    info: Info,
    samples: Samples,
}

impl Builder {
    /// Sets the chromosome.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .build();
    ///
    /// assert_eq!(record.chromosome(), "sq0");
    /// ```
    pub fn set_chromosome<C>(mut self, chromosome: C) -> Self
    where
        C: Into<String>,
    {
        self.chromosome = chromosome.into();
        self
    }

    /// Sets the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Position};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_position(Position::from(8))
    ///     .build();
    ///
    /// assert_eq!(usize::from(record.position()), 8);
    /// ```
    pub fn set_position(mut self, position: Position) -> Self {
        self.position = position;
        self
    }

    /// Sets a list of IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Ids};
    ///
    /// let ids: Ids = [String::from("nd0")].into_iter().collect();
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_ids(ids.clone())
    ///     .build();
    ///
    /// assert_eq!(record.ids(), &ids);
    /// ```
    pub fn set_ids(mut self, ids: Ids) -> Self {
        self.ids = ids;
        self
    }

    /// Sets the reference bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_reference_bases("A")
    ///     .build();
    ///
    /// assert_eq!(record.reference_bases(), "A");
    /// ```
    pub fn set_reference_bases<B>(mut self, reference_bases: B) -> Self
    where
        B: Into<String>,
    {
        self.reference_bases = reference_bases.into();
        self
    }

    /// Sets the alternate bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::AlternateBases};
    ///
    /// let alternate_bases = AlternateBases::from(vec![String::from("C")]);
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_alternate_bases(alternate_bases.clone())
    ///     .build();
    ///
    /// assert_eq!(record.alternate_bases(), &alternate_bases);
    /// ```
    pub fn set_alternate_bases(mut self, alternate_bases: AlternateBases) -> Self {
        self.alternate_bases = alternate_bases;
        self
    }

    /// Sets the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_quality_score(13.0)
    ///     .build();
    ///
    /// assert_eq!(record.quality_score(), Some(13.0));
    /// ```
    pub fn set_quality_score(mut self, quality_score: f32) -> Self {
        self.quality_score = Some(quality_score);
        self
    }

    /// Sets the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Filters};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_filters(Filters::Pass)
    ///     .build();
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// ```
    pub fn set_filters(mut self, filters: Filters) -> Self {
        self.filters = Some(filters);
        self
    }

    /// Sets additional information.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     variant::record_buf::{
    ///         info::field::{key, Value},
    ///         Info,
    ///     },
    /// };
    ///
    /// let info: Info = [
    ///     (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(3))),
    ///     (String::from(key::ALLELE_FREQUENCIES), Some(Value::from(vec![Some(0.5)]))),
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_info(info.clone())
    ///     .build();
    ///
    /// assert_eq!(record.info(), &info);
    /// ```
    pub fn set_info(mut self, info: Info) -> Self {
        self.info = info;
        self
    }

    /// Sets the list of genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     variant::record_buf::{
    ///         samples::{keys::key, sample::Value, Keys},
    ///         Samples,
    ///     },
    /// };
    ///
    /// let keys = Keys::try_from(vec![
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ])?;
    ///
    /// let samples = Samples::new(
    ///     keys,
    ///     vec![vec![Some(Value::from("0|0")), Some(Value::from(13))]],
    /// );
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_samples(samples.clone())
    ///     .build();
    ///
    /// assert_eq!(record.samples(), &samples);
    /// # Ok::<_, vcf::variant::record_buf::samples::keys::TryFromKeyVectorError>(())
    /// ```
    pub fn set_samples(mut self, samples: Samples) -> Self {
        self.samples = samples;
        self
    }

    /// Builds a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let record = vcf::variant::RecordBuf::builder().build();
    /// ```
    pub fn build(self) -> RecordBuf {
        RecordBuf {
            chromosome: self.chromosome,
            position: self.position,
            ids: self.ids,
            reference_bases: self.reference_bases,
            alternate_bases: self.alternate_bases,
            quality_score: self.quality_score,
            filters: self.filters,
            info: self.info,
            samples: self.samples,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            chromosome: String::new(),
            position: Position::from(1),
            ids: Ids::default(),
            reference_bases: String::new(),
            alternate_bases: AlternateBases::default(),
            quality_score: None,
            filters: None,
            info: Info::default(),
            samples: Samples::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let record = Builder::default();

        assert!(record.chromosome.is_empty());
        assert_eq!(record.position, Position::from(1));
        assert!(record.ids.is_empty());
        assert!(record.reference_bases.is_empty());
        assert!(record.alternate_bases.is_empty());
        assert!(record.quality_score.is_none());
        assert!(record.filters.is_none());
        assert!(record.info.is_empty());
        assert!(record.samples.is_empty());
    }
}
