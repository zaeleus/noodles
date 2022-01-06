//! VCF record and fields.

pub mod alternate_bases;
pub mod builder;
pub mod chromosome;
mod field;
pub mod filters;
pub mod genotypes;
pub mod ids;
pub mod info;
mod parser;
pub mod position;
pub mod quality_score;
pub mod reference_bases;
pub(crate) mod value;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, chromosome::Chromosome, field::Field,
    filters::Filters, genotypes::Genotypes, ids::Ids, info::Info, parser::ParseError,
    position::Position, quality_score::QualityScore, reference_bases::ReferenceBases,
};

#[deprecated(
    since = "0.8.0",
    note = "Use `noodles_vcf::record::genotypes::{genotype, Genotype}` instead."
)]
pub use self::genotypes::{genotype, Genotype};

#[deprecated(
    since = "0.11.0",
    note = "Use `noodles_vcf::record::genotypes::Format` instead."
)]
pub use self::genotypes::Keys as Format;

use std::{error, fmt, num, str::FromStr};

use super::Header;

pub(crate) const MISSING_FIELD: &str = ".";
pub(crate) const FIELD_DELIMITER: char = '\t';

/// A VCF record.
///
/// A VCF record has 8 required fields: chromosome (`CHROM`), position (`POS`), IDs (`ID`),
/// reference bases (`REF`), alternate bases (`ALT`), quality score (`QUAL`), filters (`FILTER`),
/// and information (`INFO`).
///
/// Additionally, each record can have genotype information. This adds the extra `FORMAT` field and
/// a number of genotype fields.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    chromosome: Chromosome,
    position: Position,
    ids: Ids,
    reference_bases: ReferenceBases,
    alternate_bases: AlternateBases,
    quality_score: Option<QualityScore>,
    filters: Option<Filters>,
    info: Info,
    genotypes: Genotypes,
}

impl Record {
    /// Parses a raw VCF record.
    pub fn try_from_str(s: &str, header: &Header) -> Result<Self, ParseError> {
        parser::parse(s, header)
    }

    /// Returns a builder to create a record from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let builder = vcf::Record::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the chromosome of the record.
    ///
    /// The chromosome is either a reference sequence name or a symbol (`<identifier>`).
    ///
    /// This is a required field and guaranteed to be set.
    ///
    /// # Examples
    ///
    /// ```
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
    pub fn chromosome(&self) -> &Chromosome {
        &self.chromosome
    }

    /// Returns a mutable reference to the chromosome.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Chromosome, Position}};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.chromosome_mut() = "sq1".parse()?;
    ///
    /// assert_eq!(record.chromosome().to_string(), "sq1");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn chromosome_mut(&mut self) -> &mut Chromosome {
        &mut self.chromosome
    }

    /// Returns the start position of the reference bases or indicates a telomeric breakend.
    ///
    /// This field is overloaded. If the record represents a telomere, the telomeric breakends are
    /// set to 0 and _n_ + 1, where _n_ is the length of the chromosome. Otherwise, it is a 1-based
    /// start position of the reference bases.
    ///
    /// This is a required field. It is guaranteed to be set and >= 0.
    ///
    /// # Examples
    ///
    /// ```
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
    pub fn position(&self) -> Position {
        self.position
    }

    /// Returns a mutable reference to the start position of the reference bases.
    ///
    /// This field is overloaded. If the record represents a telomere, the telomeric breakends are
    /// set to 0 and _n_ + 1, where _n_ is the length of the chromosome. Otherwise, it is a 1-based
    /// start position of the reference bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(8)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.position_mut() = Position::try_from(13)?;
    ///
    /// assert_eq!(i32::from(record.position()), 13);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn position_mut(&mut self) -> &mut Position {
        &mut self.position
    }

    /// Returns a list of IDs of the record.
    ///
    /// # Examples
    ///
    /// ```
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
    pub fn ids(&self) -> &Ids {
        &self.ids
    }

    /// Returns a mutable reference to the IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Ids, Position}};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// let ids: Ids = "nd0".parse()?;
    /// *record.ids_mut() = ids.clone();
    ///
    /// assert_eq!(record.ids(), &ids);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn ids_mut(&mut self) -> &mut Ids {
        &mut self.ids
    }

    /// Returns the reference bases of the record.
    ///
    /// This is a required field and guaranteed to be nonempty.
    ///
    /// # Examples
    ///
    /// ```
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
    pub fn reference_bases(&self) -> &ReferenceBases {
        &self.reference_bases
    }

    /// Returns a mutable reference to the reference bases of the record.
    ///
    /// This is a required field and guaranteed to be nonempty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{reference_bases::Base, Position, ReferenceBases},
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.reference_bases_mut() = "T".parse()?;
    ///
    /// assert_eq!(
    ///     record.reference_bases(),
    ///     &ReferenceBases::try_from(vec![Base::T])?,
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_bases_mut(&mut self) -> &mut ReferenceBases {
        &mut self.reference_bases
    }

    /// Returns the alternate bases of the record.
    ///
    /// # Examples
    ///
    /// ```
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
    pub fn alternate_bases(&self) -> &AlternateBases {
        &self.alternate_bases
    }

    /// Returns a mutable reference to the alternate bases of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{alternate_bases::Allele, reference_bases::Base, AlternateBases, Position},
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.alternate_bases_mut() = "C".parse()?;
    ///
    /// assert_eq!(
    ///     record.alternate_bases(),
    ///     &AlternateBases::from(vec![Allele::Bases(vec![Base::C])]),
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternate_bases_mut(&mut self) -> &mut AlternateBases {
        &mut self.alternate_bases
    }

    /// Returns the quality score of the record.
    ///
    /// The quality score is a [Phred quality score].
    ///
    /// [Phred quality score]: https://en.wikipedia.org/wiki/Phred_quality_score
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Position, QualityScore}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_quality_score(QualityScore::try_from(13.0)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.quality_score().map(f32::from), Some(13.0));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn quality_score(&self) -> Option<QualityScore> {
        self.quality_score
    }

    /// Returns a mutable reference to the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Position, QualityScore}};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.quality_score_mut() = QualityScore::try_from(13.0).map(Some)?;
    ///
    /// assert_eq!(record.quality_score().map(f32::from), Some(13.0));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn quality_score_mut(&mut self) -> &mut Option<QualityScore> {
        &mut self.quality_score
    }

    /// Returns the filters of the record.
    ///
    /// The filters can either be pass (`PASS`), a list of filter names that caused the record to
    /// fail, (e.g., `q10`), or missing (`.`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Filters, Position}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_filters(Filters::Pass)
    ///     .build()?;
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn filters(&self) -> Option<&Filters> {
        self.filters.as_ref()
    }

    /// Returns a mutable reference to the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Filters, Position}};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.filters_mut() = Some(Filters::Pass);
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn filters_mut(&mut self) -> &mut Option<Filters> {
        &mut self.filters
    }

    /// Returns the addition information of the record.
    ///
    /// # Examples
    ///
    /// ```
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
    ///     Field::new(Key::SamplesWithDataCount, Some(Value::Integer(3))),
    ///     Field::new(Key::AlleleFrequencies, Some(Value::FloatArray(vec![0.5]))),
    /// ])?;
    ///
    /// assert_eq!(record.info(), &expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn info(&self) -> &Info {
        &self.info
    }

    /// Returns a mutable reference to the additional info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{info::{field::{Key, Value}, Field}, Info, Position},
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_alternate_bases("C".parse()?)
    ///     .set_info("NS=3;AF=0.5".parse()?)
    ///     .build()?;
    ///
    /// let dp = Field::new(Key::TotalDepth, Some(Value::Integer(13)));
    /// record.info_mut().insert(dp.clone());
    ///
    /// let expected = Info::try_from(vec![
    ///     Field::new(Key::SamplesWithDataCount, Some(Value::Integer(3))),
    ///     Field::new(Key::AlleleFrequencies, Some(Value::FloatArray(vec![0.5]))),
    ///     dp,
    /// ])?;
    ///
    /// assert_eq!(record.info(), &expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn info_mut(&mut self) -> &mut Info {
        &mut self.info
    }

    /// Returns the format of the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{genotypes::{genotype::field::Key, Genotype, Keys}, Genotypes, Position},
    /// };
    ///
    /// let keys = "GT:GQ".parse()?;
    /// let genotypes = vec![Genotype::from_str_format("0|0:13", &keys)?];
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_genotypes(Genotypes::new(keys, genotypes))
    ///     .build()?;
    ///
    /// assert_eq!(record.format(), &Keys::try_from(vec![
    ///     Key::Genotype,
    ///     Key::ConditionalGenotypeQuality,
    /// ])?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn format(&self) -> &genotypes::Keys {
        self.genotypes.keys()
    }

    /// Returns the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{
    ///         genotypes::{genotype::{field::{Key, Value}, Field}, Genotype, Genotypes},
    ///         Position,
    ///     },
    /// };
    ///
    /// let keys = "GT:GQ".parse()?;
    /// let genotypes = Genotypes::new(
    ///     keys,
    ///     vec![Genotype::try_from(vec![
    ///         Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ])?],
    /// );
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_genotypes(genotypes.clone())
    ///     .build()?;
    ///
    ///
    /// assert_eq!(record.genotypes(), &genotypes);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotypes(&self) -> &Genotypes {
        &self.genotypes
    }

    /// Returns a mutable reference to the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{
    ///         genotypes::{genotype::{field::{Key, Value}, Field}, Genotype, Genotypes},
    ///         Position,
    ///     },
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// let keys = "GT:GQ".parse()?;
    /// let genotypes = Genotypes::new(
    ///     keys,
    ///     vec![Genotype::try_from(vec![
    ///         Field::new(Key::Genotype, Some(Value::String(String::from("0|0")))),
    ///         Field::new(Key::ConditionalGenotypeQuality, Some(Value::Integer(13))),
    ///     ])?],
    /// );
    ///
    /// *record.genotypes_mut() = genotypes.clone();
    ///
    /// assert_eq!(record.genotypes(), &genotypes);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotypes_mut(&mut self) -> &mut Genotypes {
        &mut self.genotypes
    }
}

/// An error returned when the end position is invalid.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EndError {
    /// The position is invalid.
    InvalidPosition(position::TryFromIntError),
    /// The INFO end position (`END`) field value type is invalid.
    InvalidInfoEndPositionFieldValue,
    /// The reference bases length is invalid (> [`i32::MAX`]).
    InvalidReferenceBasesLength(num::TryFromIntError),
    /// The calculation of the end position overflowed.
    PositionOverflow(i32, i32),
}

impl error::Error for EndError {}

impl fmt::Display for EndError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPosition(e) => write!(f, "invalid position: {}", e),
            Self::InvalidInfoEndPositionFieldValue => {
                write!(f, "invalid INFO end position (`END`) field value type")
            }
            Self::InvalidReferenceBasesLength(e) => {
                write!(f, "invalid reference base length: {}", e)
            }
            Self::PositionOverflow(start, len) => write!(
                f,
                "calculation of the end position overflowed: {} + {}",
                start, len,
            ),
        }
    }
}

impl Record {
    /// Returns or calculates the end position on the reference sequence.
    ///
    /// If available, this returns the value of the `END` INFO field. Otherwise, it is calculated
    /// using the start position and reference bases length.
    ///
    /// The end position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ## From the `END` INFO field value
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Chromosome, Position}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("ACGT".parse()?)
    ///     .set_info("END=8".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.end(), Ok(Position::try_from(8)?));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// ## Calculated using the start position and reference bases length
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Chromosome, Position}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("ACGT".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.end(), Ok(Position::try_from(4)?));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> Result<Position, EndError> {
        use info::field::{Key, Value};

        let end = if let Some(field) = self.info().get(&Key::EndPosition) {
            match field.value() {
                Some(Value::Integer(n)) => *n,
                _ => return Err(EndError::InvalidInfoEndPositionFieldValue),
            }
        } else {
            let start = i32::from(self.position());
            // `len` is guaranteed to be > 0.
            let len = i32::try_from(self.reference_bases().len())
                .map_err(EndError::InvalidReferenceBasesLength)?;
            start
                .checked_add(len - 1)
                .ok_or(EndError::PositionOverflow(start, len))?
        };

        Position::try_from(end).map_err(EndError::InvalidPosition)
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{chrom}\t{pos}\t{id}\t{ref}\t{alt}",
            chrom = self.chromosome(),
            pos = i32::from(self.position()),
            id = self.ids(),
            r#ref = self.reference_bases(),
            alt = self.alternate_bases(),
        )?;

        if let Some(quality_score) = self.quality_score() {
            write!(f, "\t{}", quality_score)?;
        } else {
            write!(f, "\t{}", MISSING_FIELD)?;
        }

        if let Some(filters) = self.filters() {
            write!(f, "\t{}", filters)?;
        } else {
            write!(f, "\t{}", MISSING_FIELD)?;
        }

        write!(f, "\t{}", self.info())?;

        if !self.genotypes().is_empty() {
            write!(f, "\t{}", self.genotypes())?;
        }

        Ok(())
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from_str(s, &Header::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_end() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .set_info(Info::try_from(vec![info::Field::new(
                info::field::Key::EndPosition,
                Some(info::field::Value::Flag),
            )])?)
            .build()?;

        assert_eq!(
            record.end(),
            Err(EndError::InvalidInfoEndPositionFieldValue)
        );

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(i32::MAX)?)
            .set_reference_bases("ACGT".parse()?)
            .build()?;

        assert_eq!(record.end(), Err(EndError::PositionOverflow(i32::MAX, 4)));

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .build()?;

        assert_eq!(record.to_string(), "sq0\t1\t.\tA\t.\t.\t.\t.");

        Ok(())
    }

    #[test]
    fn test_fmt_with_format() -> Result<(), Box<dyn std::error::Error>> {
        use genotypes::Genotype;

        let keys: genotypes::Keys = "GT:GQ".parse()?;
        let genotypes = vec![Genotype::from_str_format("0|0:13", &keys)?];

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .set_genotypes(Genotypes::new(keys, genotypes))
            .build()?;

        assert_eq!(
            record.to_string(),
            "sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13"
        );

        Ok(())
    }
}
