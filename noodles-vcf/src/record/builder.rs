//! VCF record builder.

use std::{convert::TryFrom, error, fmt};

use super::{
    reference_bases::Base, AlternateBases, Chromosome, Filters, Format, Genotype, Genotypes, Ids,
    Info, Position, QualityScore, Record, ReferenceBases,
};

/// A VCF record builder.
#[derive(Debug, Default, PartialEq)]
pub struct Builder {
    chromosome: Option<Chromosome>,
    position: Option<Position>,
    ids: Ids,
    reference_bases: Vec<Base>,
    alternate_bases: AlternateBases,
    quality_score: QualityScore,
    filters: Filters,
    info: Info,
    format: Option<Format>,
    genotypes: Genotypes,
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
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::{Chromosome, Position}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.chromosome(), &Chromosome::Name(String::from("sq0")));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_chromosome(mut self, chromosome: Chromosome) -> Self {
        self.chromosome = Some(chromosome);
        self
    }

    /// Sets the start position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(8)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(i32::from(record.position()), 8);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_ids("nd0".parse()?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(*record.ids(), "nd0".parse()?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{reference_bases::Base, Position, ReferenceBases},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     record.reference_bases(),
    ///     &ReferenceBases::try_from(vec![Base::A])?,
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_reference_bases(mut self, reference_bases: ReferenceBases) -> Self {
        self.reference_bases.clear();
        self.reference_bases.extend(reference_bases.iter());
        self
    }

    /// Adds a base to reference bases.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{reference_bases::Base, Position, ReferenceBases},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .add_reference_base(Base::A)
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     record.reference_bases(),
    ///     &ReferenceBases::try_from(vec![Base::A])?,
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_reference_base(mut self, reference_base: Base) -> Self {
        self.reference_bases.push(reference_base);
        self
    }

    /// Sets the alternate bases.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{alternate_bases::Allele, reference_bases::Base, AlternateBases, Position},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_alternate_bases("C".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     record.alternate_bases(),
    ///     &AlternateBases::from(vec![Allele::Bases(vec![Base::C])]),
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::{Position, QualityScore}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_quality_score(QualityScore::try_from(13.0)?)
    ///     .build()?;
    ///
    /// assert_eq!(*record.quality_score(), Some(13.0));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_quality_score(mut self, quality_score: QualityScore) -> Self {
        self.quality_score = quality_score;
        self
    }

    /// Sets the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::{Filters, Position}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_filters(Filters::Pass)
    ///     .build()?;
    ///
    /// assert_eq!(record.filters(), &Filters::Pass);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_filters(mut self, filters: Filters) -> Self {
        self.filters = filters;
        self
    }

    /// Sets additional information.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{info::{field::{Key, Value}, Field}, Info, Position},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_alternate_bases("C".parse()?)
    ///     .set_info("NS=3;AF=0.5".parse()?)
    ///     .build()?;
    ///
    /// let expected = Info::try_from(vec![
    ///     Field::new(Key::SamplesWithDataCount, Value::Integer(3)),
    ///     Field::new(Key::AlleleFrequencies, Value::FloatArray(vec![0.5])),
    /// ])?;
    ///
    /// assert_eq!(record.info(), &expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_info(mut self, info: Info) -> Self {
        self.info = info;
        self
    }

    /// Sets the genotype format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{genotype::field::Key, Format, Genotype, Position},
    /// };
    ///
    /// let format: Format = "GT:GQ".parse()?;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_format(format.clone())
    ///     .add_genotype(Genotype::from_str_format("0|0:13", &format)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.format(), Some(&Format::try_from(vec![
    ///     Key::Genotype,
    ///     Key::ConditionalGenotypeQuality,
    /// ])?));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the list of genotypes.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{genotype::{field::{Key, Value}, Field, Genotypes}, Format, Genotype, Position},
    /// };
    ///
    /// let format: Format = "GT:GQ".parse()?;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_format(format.clone())
    ///     .set_genotypes(vec![Genotype::from_str_format("0|0:13", &format)?].into())
    ///     .build()?;
    ///
    /// let expected = vec![
    ///     Genotype::try_from(vec![
    ///         Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ])?,
    /// ].into();
    ///
    /// assert_eq!(record.genotypes(), &expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_genotypes(mut self, genotypes: Genotypes) -> Self {
        self.genotypes = genotypes;
        self
    }

    /// Adds a genotype.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{genotype::{field::{Key, Value}, Field}, Format, Genotype, Position},
    /// };
    ///
    /// let format: Format = "GT:GQ".parse()?;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_format(format.clone())
    ///     .add_genotype(Genotype::from_str_format("0|0:13", &format)?)
    ///     .build()?;
    ///
    /// let expected = vec![
    ///     Genotype::try_from(vec![
    ///         Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ])?,
    /// ].into();
    ///
    /// assert_eq!(record.genotypes(), &expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_genotype(mut self, genotype: Genotype) -> Self {
        self.genotypes.push(genotype);
        self
    }

    /// Builds a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let record = vcf::Record::builder().build();
    /// ```
    pub fn build(self) -> Result<Record, BuildError> {
        Ok(Record {
            chromosome: self.chromosome.ok_or(BuildError::MissingChromosome)?,
            position: self.position.ok_or(BuildError::MissingPosition)?,
            ids: self.ids,
            reference_bases: ReferenceBases::try_from(self.reference_bases)
                .map_err(|_| BuildError::MissingReferenceBases)?,
            alternate_bases: self.alternate_bases,
            quality_score: self.quality_score,
            filters: self.filters,
            info: self.info,
            format: self.format,
            genotypes: self.genotypes,
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
        assert!(matches!(record.filters, Filters::Missing));
        assert!(record.info.is_empty());
        assert!(record.format.is_none());
        assert!(record.genotypes.is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let result = Builder::default()
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .build();
        assert_eq!(result, Err(BuildError::MissingChromosome));

        let result = Builder::default()
            .set_chromosome("sq0".parse()?)
            .set_reference_bases("A".parse()?)
            .build();
        assert_eq!(result, Err(BuildError::MissingPosition));

        let result = Builder::default()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .build();
        assert_eq!(result, Err(BuildError::MissingReferenceBases));

        Ok(())
    }
}
