//! VCF record builder.

use std::{error, fmt};

use super::{AlternateBases, Filters, Ids, Info, Position, RecordBuf, Samples};

/// A VCF record builder.
#[derive(Debug, Default, PartialEq)]
pub struct Builder {
    chromosome: Option<String>,
    position: Option<Position>,
    ids: Ids,
    reference_bases: String,
    alternate_bases: AlternateBases,
    quality_score: Option<f32>,
    filters: Option<Filters>,
    info: Info,
    samples: Samples,
}

/// An error returned when a VCF record fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The chromosome is missing.
    MissingChromosome,
    /// The position is missing.
    MissingPosition,
    /// The reference bases are missing.
    MissingReferenceBases,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingChromosome => f.write_str("missing chromosome"),
            Self::MissingPosition => f.write_str("missing position"),
            Self::MissingReferenceBases => f.write_str("missing reference bases"),
        }
    }
}

impl Builder {
    /// Sets the chromosome.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Position};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(record.chromosome(), "sq0");
    /// # Ok::<_, vcf::variant::record_buf::builder::BuildError>(())
    /// ```
    pub fn set_chromosome<C>(mut self, chromosome: C) -> Self
    where
        C: Into<String>,
    {
        self.chromosome = Some(chromosome.into());
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(8))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(usize::from(record.position()), 8);
    /// # Ok::<_, vcf::variant::record_buf::builder::BuildError>(())
    /// ```
    pub fn set_position(mut self, position: Position) -> Self {
        self.position = Some(position);
        self
    }

    /// Sets a list of IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::{Ids, Position}};
    ///
    /// let ids: Ids = [String::from("nd0")].into_iter().collect();
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_ids(ids.clone())
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(record.ids(), &ids);
    /// # Ok::<(), vcf::variant::record_buf::builder::BuildError>(())
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
    /// use noodles_vcf::{self as vcf, variant::record_buf::Position};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(record.reference_bases(), "A");
    /// # Ok::<_, vcf::variant::record_buf::builder::BuildError>(())
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
    /// use noodles_vcf::{self as vcf, variant::record_buf::{AlternateBases, Position}};
    ///
    /// let alternate_bases = AlternateBases::from(vec![String::from("C")]);
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_alternate_bases(alternate_bases.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.alternate_bases(), &alternate_bases);
    /// # Ok::<_, vcf::variant::record_buf::builder::BuildError>(())
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
    /// use noodles_vcf::{self as vcf, variant::record_buf::Position};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_quality_score(13.0)
    ///     .build()?;
    ///
    /// assert_eq!(record.quality_score(), Some(13.0));
    /// # Ok::<_, vcf::variant::record_buf::builder::BuildError>(())
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
    /// use noodles_vcf::{self as vcf, variant::record_buf::{Filters, Position}};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_filters(Filters::Pass)
    ///     .build()?;
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// # Ok::<_, vcf::variant::record_buf::builder::BuildError>(())
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
    ///     variant::record_buf::{info::field::{key, Value}, Info, Position},
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_info(info.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.info(), &info);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    ///         Position, Samples,
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_samples(samples.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.samples(), &samples);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    pub fn build(self) -> Result<RecordBuf, BuildError> {
        Ok(RecordBuf {
            chromosome: self.chromosome.ok_or(BuildError::MissingChromosome)?,
            position: self.position.ok_or(BuildError::MissingPosition)?,
            ids: self.ids,
            reference_bases: if self.reference_bases.is_empty() {
                return Err(BuildError::MissingReferenceBases);
            } else {
                self.reference_bases
            },
            alternate_bases: self.alternate_bases,
            quality_score: self.quality_score,
            filters: self.filters,
            info: self.info,
            samples: self.samples,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let record = Builder::default();

        assert!(record.chromosome.is_none());
        assert!(record.position.is_none());
        assert!(record.ids.is_empty());
        assert!(record.reference_bases.is_empty());
        assert!(record.alternate_bases.is_empty());
        assert!(record.quality_score.is_none());
        assert!(record.filters.is_none());
        assert!(record.info.is_empty());
        assert!(record.samples.is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let result = Builder::default()
            .set_position(Position::from(1))
            .set_reference_bases("A")
            .build();
        assert_eq!(result, Err(BuildError::MissingChromosome));

        let result = Builder::default()
            .set_chromosome("sq0")
            .set_reference_bases("A")
            .build();
        assert_eq!(result, Err(BuildError::MissingPosition));

        let result = Builder::default()
            .set_chromosome("sq0")
            .set_position(Position::from(1))
            .build();
        assert_eq!(result, Err(BuildError::MissingReferenceBases));

        Ok(())
    }
}
