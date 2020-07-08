//! VCF record and fields.

pub mod alternate_bases;
mod builder;
mod chromosome;
mod field;
mod filter_status;
pub mod format;
pub mod genotype;
mod id;
pub mod info;
mod quality_score;
pub mod reference_bases;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, chromosome::Chromosome, field::Field,
    filter_status::FilterStatus, format::Format, genotype::Genotype, id::Id, info::Info,
    quality_score::QualityScore, reference_bases::ReferenceBases,
};

use std::{error, fmt, num, str::FromStr};

pub(crate) const MISSING_FIELD: &str = ".";
const FIELD_DELIMITER: char = '\t';

/// A VCF record.
///
/// A VCF record has 8 required fields:
///
///   1. chromosome (`CRHOM`),
///   2. position (`POS`),
///   3. id (`ID`),
///   4. reference bases (`REF`),
///   5. alternate bases (`ALT`),
///   6. quality score (`QUAL`),
///   7. filter status (`FILTER`), and
///   8. information (`INFO`),
///
/// Additionally, each record can have genotype information. This adds the extra `FORMAT` field
/// and a number of genotype fields.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    chromosome: Chromosome,
    position: i32,
    id: Id,
    reference_bases: ReferenceBases,
    alternate_bases: AlternateBases,
    quality_score: QualityScore,
    filter_status: FilterStatus,
    info: Info,
    format: Option<Format>,
    genotypes: Vec<Genotype>,
}

impl Record {
    /// Returns a builder to create a record from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let builder = vcf::Record::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::new()
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
    /// use noodles_vcf::{self as vcf, record::Chromosome};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(1)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.chromosome(), &Chromosome::Name(String::from("sq0")));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn chromosome(&self) -> &Chromosome {
        &self.chromosome
    }

    /// Returns the start position of the record.
    ///
    /// VCF positions are 1-based.
    ///
    /// This is a required field and guaranteed to be set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(8)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.position(), 8);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn position(&self) -> i32 {
        self.position
    }

    /// Returns a list of IDs of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(8)
    ///     .set_id("nd0".parse()?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(**record.id(), [String::from("nd0")]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn id(&self) -> &Id {
        &self.id
    }

    /// Returns the reference bases of the record.
    ///
    /// This is a required field and guaranteed to be nonempty.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::{reference_bases::Base, ReferenceBases}};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(8)
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

    /// Returns the alternate bases of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{alternate_bases::Allele, reference_bases::Base, AlternateBases}
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(8)
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

    /// Returns the quality score of the record.
    ///
    /// The quality score is a [Phred quality score].
    ///
    /// [Phred quality score]: https://en.wikipedia.org/wiki/Phred_quality_score
    ///
    /// # Examples
    ///
    /// ```
    /// use std::convert::TryFrom;
    /// use noodles_vcf::{self as vcf, record::QualityScore};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(8)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_quality_score(QualityScore::try_from(13.0)?)
    ///     .build()?;
    ///
    /// assert_eq!(*record.quality_score(), Some(13.0));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn quality_score(&self) -> QualityScore {
        self.quality_score
    }

    /// Returns the filter status of the record.
    ///
    /// The filter status can either be pass (`PASS`), a list of filter names
    /// that caused the record to fail, (e.g., `q10`), or missing (`.`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::FilterStatus};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(1)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_filter_status(FilterStatus::Pass)
    ///     .build()?;
    ///
    /// assert_eq!(record.filter_status(), &FilterStatus::Pass);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn filter_status(&self) -> &FilterStatus {
        &self.filter_status
    }

    /// Returns the addition information of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{info::{field::{Key, Value}, Field}, Info},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(1)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_alternate_bases("C".parse()?)
    ///     .set_info("NS=3;AF=0.5".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.info(), &Info::from(vec![
    ///     Field::new(Key::SamplesWithDataCount, Value::Integer(3)),
    ///     Field::new(Key::AlleleFrequencies, Value::FloatArray(vec![0.5])),
    /// ]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn info(&self) -> &Info {
        &self.info
    }

    /// Returns the format of the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{genotype::field::Key, Format, Genotype}
    /// };
    ///
    /// let format: Format = "GT:GQ".parse()?;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(1)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_format(format.clone())
    ///     .add_genotype(Genotype::from_str_format("0|0:13", &format)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.format(), Some(&Format::from(vec![
    ///     Key::Genotype,
    ///     Key::ConditionalGenotypeQuality,
    /// ])));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn format(&self) -> Option<&Format> {
        self.format.as_ref()
    }

    /// Returns the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{genotype::{field::{Key, Value}, Field}, Format, Genotype}
    /// };
    ///
    /// let format: Format = "GT:GQ".parse()?;
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(1)
    ///     .set_reference_bases("A".parse()?)
    ///     .set_format(format.clone())
    ///     .add_genotype(Genotype::from_str_format("0|0:13", &format)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.genotypes(), [Genotype::from(vec![
    ///     Field::new(Key::Genotype, Value::String(String::from("0|0"))),
    ///     Field::new(Key::ConditionalGenotypeQuality, Value::Integer(13)),
    /// ])]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotypes(&self) -> &[Genotype] {
        &self.genotypes
    }
}

/// An error returned when a raw VCF record fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Field),
    /// The chromosome is invalid.
    InvalidChromosome(chromosome::ParseError),
    /// The position is invalid.
    InvalidPosition(num::ParseIntError),
    /// The ID is invalid.
    InvalidId(id::ParseError),
    /// The reference bases are invalid.
    InvalidReferenceBases(reference_bases::ParseError),
    /// The alternate bases are invalid.
    InvalidAlternateBases(alternate_bases::ParseError),
    /// The quality score is invalid.
    InvalidQualityScore(quality_score::ParseError),
    /// The filter status is invalid.
    InvalidFilterStatus(filter_status::ParseError),
    /// The info is invalid.
    InvalidInfo(info::ParseError),
    /// The format is invalid.
    InvalidFormat(format::ParseError),
    /// A genotype is invalid.
    InvalidGenotype(genotype::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(field) => write!(f, "missing field: {}", field),
            Self::InvalidChromosome(e) => write!(f, "{}", e),
            Self::InvalidPosition(e) => write!(f, "{}", e),
            Self::InvalidId(e) => write!(f, "{}", e),
            Self::InvalidReferenceBases(e) => write!(f, "{}", e),
            Self::InvalidAlternateBases(e) => write!(f, "{}", e),
            Self::InvalidQualityScore(e) => write!(f, "{}", e),
            Self::InvalidFilterStatus(e) => write!(f, "{}", e),
            Self::InvalidInfo(e) => write!(f, "{}", e),
            Self::InvalidFormat(e) => write!(f, "{}", e),
            Self::InvalidGenotype(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(FIELD_DELIMITER);

        let chrom = parse_string(&mut fields, Field::Chromosome)
            .and_then(|s| s.parse().map_err(ParseError::InvalidChromosome))?;

        let pos = parse_string(&mut fields, Field::Position)
            .and_then(|s| s.parse().map_err(ParseError::InvalidPosition))?;

        let id = parse_string(&mut fields, Field::Id)
            .and_then(|s| s.parse().map_err(ParseError::InvalidId))?;

        let r#ref = parse_string(&mut fields, Field::ReferenceBases)
            .and_then(|s| s.parse().map_err(ParseError::InvalidReferenceBases))?;

        let alt = parse_string(&mut fields, Field::AlternateBases)
            .and_then(|s| s.parse().map_err(ParseError::InvalidAlternateBases))?;

        let qual = parse_string(&mut fields, Field::QualityScore)
            .and_then(|s| s.parse().map_err(ParseError::InvalidQualityScore))?;

        let filter = parse_string(&mut fields, Field::FilterStatus)
            .and_then(|s| s.parse().map_err(ParseError::InvalidFilterStatus))?;

        let info = parse_string(&mut fields, Field::Info)
            .and_then(|s| s.parse().map_err(ParseError::InvalidInfo))?;

        let format = match fields.next() {
            Some(s) => s.parse().map(Some).map_err(ParseError::InvalidFormat)?,
            None => None,
        };

        let genotypes = format
            .as_ref()
            .map(|f| {
                fields
                    .map(|s| Genotype::from_str_format(s, f))
                    .collect::<Result<_, _>>()
                    .map_err(ParseError::InvalidGenotype)
            })
            .unwrap_or_else(|| Ok(Vec::new()))?;

        Ok(Self {
            chromosome: chrom,
            position: pos,
            id,
            reference_bases: r#ref,
            alternate_bases: alt,
            quality_score: qual,
            filter_status: filter,
            info,
            format,
            genotypes,
        })
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<&'a str, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().ok_or_else(|| ParseError::MissingField(field))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        use alternate_bases::Allele;
        use reference_bases::Base;

        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL";
        let record: Record = s.parse()?;

        assert!(matches!(record.chromosome(), Chromosome::Name(name) if name == "chr1"));

        assert_eq!(record.position(), 13);
        assert_eq!(**record.id(), [String::from("nd0")]);

        let reference_bases = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(&record.reference_bases()[..], &reference_bases[..]);

        let alternate_bases = [Allele::Bases(vec![Base::A])];
        assert_eq!(&record.alternate_bases()[..], &alternate_bases[..]);

        assert_eq!(*record.quality_score(), Some(5.8));
        assert_eq!(record.filter_status(), &FilterStatus::Pass);
        assert_eq!(record.info().len(), 1);
        assert!(record.format().is_none());
        assert!(record.genotypes().is_empty());

        Ok(())
    }

    #[test]
    fn test_from_str_with_genotype_info() -> Result<(), ParseError> {
        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL\tGT:GQ\t0|1:13";
        let record: Record = s.parse()?;

        let expected = [
            genotype::field::Key::Genotype,
            genotype::field::Key::ConditionalGenotypeQuality,
        ];

        assert_eq!(record.format().map(|f| &f[..]), Some(&expected[..]));

        let genotypes = record.genotypes();
        let expected = vec![
            genotype::Field::new(
                genotype::field::Key::Genotype,
                genotype::field::Value::String(String::from("0|1")),
            ),
            genotype::Field::new(
                genotype::field::Key::ConditionalGenotypeQuality,
                genotype::field::Value::Integer(13),
            ),
        ];

        assert_eq!(genotypes.len(), 1);
        assert_eq!(&genotypes[0][..], &expected[..]);

        Ok(())
    }
}
