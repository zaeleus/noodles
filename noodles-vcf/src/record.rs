//! VCF record and fields.

pub mod alternate_bases;
pub mod builder;
pub mod chromosome;
mod field;
pub mod filters;
pub mod format;
pub mod genotypes;
pub mod ids;
pub mod info;
pub mod position;
pub mod quality_score;
pub mod reference_bases;
pub(crate) mod value;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, chromosome::Chromosome, field::Field,
    filters::Filters, format::Format, genotypes::Genotypes, ids::Ids, info::Info,
    position::Position, quality_score::QualityScore, reference_bases::ReferenceBases,
};

#[deprecated(
    since = "0.8.0",
    note = "Use `noodles_vcf::record::genotypes::{genotype, Genotype}` instead."
)]
pub use self::genotypes::{genotype, Genotype};

use std::{convert::TryFrom, error, fmt, num, str::FromStr};

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
    quality_score: QualityScore,
    filters: Filters,
    info: Info,
    format: Option<Format>,
    genotypes: Genotypes,
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
    pub fn chromosome(&self) -> &Chromosome {
        &self.chromosome
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
    pub fn position(&self) -> Position {
        self.position
    }

    /// Returns a list of IDs of the record.
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
    pub fn ids(&self) -> &Ids {
        &self.ids
    }

    /// Returns the reference bases of the record.
    ///
    /// This is a required field and guaranteed to be nonempty.
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
    pub fn reference_bases(&self) -> &ReferenceBases {
        &self.reference_bases
    }

    /// Returns the alternate bases of the record.
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
    pub fn quality_score(&self) -> QualityScore {
        self.quality_score
    }

    /// Returns the filters of the record.
    ///
    /// The filters can either be pass (`PASS`), a list of filter names that caused the record to
    /// fail, (e.g., `q10`), or missing (`.`).
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
    pub fn filters(&self) -> &Filters {
        &self.filters
    }

    /// Returns the addition information of the record.
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
    pub fn info(&self) -> &Info {
        &self.info
    }

    /// Returns the format of the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{
    ///         genotypes::{genotype::field::Key, Genotype},
    ///         Format,
    ///         Position,
    ///     },
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
    pub fn format(&self) -> Option<&Format> {
        self.format.as_ref()
    }

    /// Returns the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{
    ///         genotypes::{genotype::{field::{Key, Value}, Field}, Genotype},
    ///         Format,
    ///         Position,
    ///     },
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
    ///     ])?
    /// ].into();
    ///
    /// assert_eq!(record.genotypes(), &expected);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn genotypes(&self) -> &Genotypes {
        &self.genotypes
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
    /// # use std::convert::TryFrom;
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
    /// # use std::convert::TryFrom;
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
                Value::Integer(n) => *n,
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
            "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}",
            chrom = self.chromosome(),
            pos = i32::from(self.position()),
            id = self.ids(),
            r#ref = self.reference_bases(),
            alt = self.alternate_bases(),
            qual = self.quality_score(),
            filter = self.filters(),
            info = self.info(),
        )?;

        if let Some(format) = self.format() {
            write!(f, "\t{}", format)?;
            write!(f, "\t{}", self.genotypes())?;
        }

        Ok(())
    }
}

/// An error returned when a raw VCF record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Field),
    /// The chromosome is invalid.
    InvalidChromosome(chromosome::ParseError),
    /// The position is invalid.
    InvalidPosition(position::ParseError),
    /// The ID is invalid.
    InvalidIds(ids::ParseError),
    /// The reference bases are invalid.
    InvalidReferenceBases(reference_bases::ParseError),
    /// The alternate bases are invalid.
    InvalidAlternateBases(alternate_bases::ParseError),
    /// The quality score is invalid.
    InvalidQualityScore(quality_score::ParseError),
    /// A filter is invalid.
    InvalidFilters(filters::ParseError),
    /// The info is invalid.
    InvalidInfo(info::ParseError),
    /// The format is invalid.
    InvalidFormat(format::ParseError),
    /// A genotype is invalid.
    InvalidGenotype(genotypes::genotype::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(field) => write!(f, "missing field: {}", field),
            Self::InvalidChromosome(e) => write!(f, "invalid chromosome: {}", e),
            Self::InvalidPosition(e) => write!(f, "invalid position: {}", e),
            Self::InvalidIds(e) => write!(f, "invalid IDs: {}", e),
            Self::InvalidReferenceBases(e) => write!(f, "invalid reference bases: {}", e),
            Self::InvalidAlternateBases(e) => write!(f, "invalid alternate bases: {}", e),
            Self::InvalidQualityScore(e) => write!(f, "invalid quality score: {}", e),
            Self::InvalidFilters(e) => write!(f, "invalid filters: {}", e),
            Self::InvalidInfo(e) => write!(f, "invalid info: {}", e),
            Self::InvalidFormat(e) => write!(f, "invalid format: {}", e),
            Self::InvalidGenotype(e) => write!(f, "invalid genotype: {}", e),
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

        let ids = parse_string(&mut fields, Field::Ids)
            .and_then(|s| s.parse().map_err(ParseError::InvalidIds))?;

        let r#ref = parse_string(&mut fields, Field::ReferenceBases)
            .and_then(|s| s.parse().map_err(ParseError::InvalidReferenceBases))?;

        let alt = parse_string(&mut fields, Field::AlternateBases)
            .and_then(|s| s.parse().map_err(ParseError::InvalidAlternateBases))?;

        let qual = parse_string(&mut fields, Field::QualityScore)
            .and_then(|s| s.parse().map_err(ParseError::InvalidQualityScore))?;

        let filter = parse_string(&mut fields, Field::Filters)
            .and_then(|s| s.parse().map_err(ParseError::InvalidFilters))?;

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
                    .map(|s| genotypes::Genotype::from_str_format(s, f))
                    .collect::<Result<_, _>>()
                    .map_err(ParseError::InvalidGenotype)
            })
            .unwrap_or_else(|| Ok(Vec::new()))?
            .into();

        Ok(Self {
            chromosome: chrom,
            position: pos,
            ids,
            reference_bases: r#ref,
            alternate_bases: alt,
            quality_score: qual,
            filters: filter,
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
    fields.next().ok_or(ParseError::MissingField(field))
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use super::*;

    #[test]
    fn test_end() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .set_info(Info::try_from(vec![info::Field::new(
                info::field::Key::EndPosition,
                info::field::Value::Flag,
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

        let format: Format = "GT:GQ".parse()?;

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .set_format(format.clone())
            .add_genotype(Genotype::from_str_format("0|0:13", &format)?)
            .build()?;

        assert_eq!(
            record.to_string(),
            "sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13"
        );

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        use alternate_bases::Allele;
        use reference_bases::Base;

        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL";
        let record: Record = s.parse()?;

        assert!(matches!(record.chromosome(), Chromosome::Name(name) if name == "chr1"));

        assert_eq!(i32::from(record.position()), 13);

        let ids = record.ids();
        assert_eq!(ids.len(), 1);
        assert!(ids.contains("nd0"));

        let reference_bases = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(&record.reference_bases()[..], &reference_bases[..]);

        let alternate_bases = [Allele::Bases(vec![Base::A])];
        assert_eq!(&record.alternate_bases()[..], &alternate_bases[..]);

        assert_eq!(*record.quality_score(), Some(5.8));
        assert_eq!(record.filters(), &Filters::Pass);
        assert_eq!(record.info().len(), 1);
        assert!(record.format().is_none());
        assert!(record.genotypes().is_empty());

        Ok(())
    }

    #[test]
    fn test_from_str_with_genotype_info() -> Result<(), Box<dyn std::error::Error>> {
        use genotypes::genotype;

        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL\tGT:GQ\t0|1:13";
        let record: Record = s.parse()?;

        let expected = Format::try_from(vec![
            genotype::field::Key::Genotype,
            genotype::field::Key::ConditionalGenotypeQuality,
        ])?;

        assert_eq!(record.format(), Some(&expected));

        let genotypes = record.genotypes();

        assert_eq!(genotypes.len(), 1);

        let actual: Vec<_> = genotypes[0].values().cloned().collect();
        let expected = vec![
            genotype::Field::new(
                genotype::field::Key::Genotype,
                Some(genotype::field::Value::String(String::from("0|1"))),
            ),
            genotype::Field::new(
                genotype::field::Key::ConditionalGenotypeQuality,
                Some(genotype::field::Value::Integer(13)),
            ),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
