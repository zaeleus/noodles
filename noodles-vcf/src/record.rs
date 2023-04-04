//! VCF record and fields.

pub mod alternate_bases;
pub mod builder;
pub mod chromosome;
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
    alternate_bases::AlternateBases, builder::Builder, chromosome::Chromosome, filters::Filters,
    genotypes::Genotypes, ids::Ids, info::Info, position::Position, quality_score::QualityScore,
    reference_bases::ReferenceBases,
};

use std::{error, fmt, num, str::FromStr};

use super::{reader::record::ParseError, Header};

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
    /// Parses a raw VCF record using a VCF header as context.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let s = "sq0\t8\t.\tA\t.\t.\tPASS\t.";
    /// let header = vcf::Header::default();
    /// let record = vcf::Record::try_from_str(s, &header)?;
    /// # Ok::<_, vcf::reader::record::ParseError>(())
    /// ```
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(8))
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(usize::from(record.position()), 8);
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
    ///     .set_position(Position::from(8))
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// *record.position_mut() = Position::from(13);
    ///
    /// assert_eq!(usize::from(record.position()), 13);
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     .set_position(Position::from(1))
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
    ///     header::info::key,
    ///     record::{info::field::Value, Info, Position},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .set_alternate_bases("C".parse()?)
    ///     .set_info("NS=3;AF=0.5".parse()?)
    ///     .build()?;
    ///
    /// let expected = [
    ///     (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(3))),
    ///     (key::ALLELE_FREQUENCIES, Some(Value::FloatArray(vec![Some(0.5)]))),
    /// ].into_iter().collect();
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
    ///     header::info::key,
    ///     record::{info::field::Value, Info, Position},
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .set_alternate_bases("C".parse()?)
    ///     .set_info("NS=3;AF=0.5".parse()?)
    ///     .build()?;
    ///
    /// record.info_mut().insert(key::TOTAL_DEPTH, Some(Value::Integer(13)));
    ///
    /// let expected = [
    ///     (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(3))),
    ///     (key::ALLELE_FREQUENCIES, Some(Value::FloatArray(vec![Some(0.5)]))),
    ///     (key::TOTAL_DEPTH, Some(Value::Integer(13))),
    /// ].into_iter().collect();
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
    ///     header::{format::key, record::value::{map::Format, Map}},
    ///     record::{genotypes::Keys, Genotypes, Position},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
    ///     .add_format(key::CONDITIONAL_GENOTYPE_QUALITY, Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY))
    ///     .build();
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .set_genotypes(Genotypes::parse("GT:GQ\t0|0:13", &header)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.format(), &Keys::try_from(vec![
    ///     key::GENOTYPE,
    ///     key::CONDITIONAL_GENOTYPE_QUALITY,
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
    ///     header::format::key,
    ///     record::{
    ///         genotypes::{sample::Value, Genotypes},
    ///         Position,
    ///     },
    /// };
    ///
    /// let keys = "GT:GQ".parse()?;
    /// let genotypes = Genotypes::new(
    ///     keys,
    ///     vec![vec![
    ///         Some(Value::String(String::from("0|0"))),
    ///         Some(Value::Integer(13)),
    ///     ]],
    /// );
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .set_genotypes(genotypes.clone())
    ///     .build()?;
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
    ///     header::format::key,
    ///     record::{
    ///         genotypes::{sample::Value, Genotypes},
    ///         Position,
    ///     },
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// let keys = "GT:GQ".parse()?;
    /// let genotypes = Genotypes::new(
    ///     keys,
    ///     vec![vec![
    ///         Some(Value::String(String::from("0|0"))),
    ///         Some(Value::Integer(13)),
    ///     ]],
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

impl Default for Record {
    fn default() -> Self {
        use self::reference_bases::Base;

        Self {
            chromosome: Chromosome::Name(String::from(".")),
            position: Position::from(1),
            ids: Ids::default(),
            // SAFETY: `[N]` is a valid list of reference bases.
            reference_bases: ReferenceBases::try_from(vec![Base::N]).unwrap(),
            alternate_bases: AlternateBases::default(),
            quality_score: None,
            filters: None,
            info: Info::default(),
            genotypes: Genotypes::default(),
        }
    }
}

/// An error returned when the end position is invalid.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EndError {
    /// The position is invalid.
    InvalidPosition(num::TryFromIntError),
    /// The INFO end position (`END`) field value type is invalid.
    InvalidInfoEndPositionFieldValue,
    /// The reference bases length is invalid (> [`i32::MAX`]).
    InvalidReferenceBasesLength(num::TryFromIntError),
    /// The calculation of the end position overflowed.
    PositionOverflow(usize, usize),
}

impl error::Error for EndError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidPosition(e) | Self::InvalidReferenceBasesLength(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for EndError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPosition(_) => f.write_str("invalid position"),
            Self::InvalidInfoEndPositionFieldValue => {
                write!(f, "invalid INFO end position (`END`) field value type")
            }
            Self::InvalidReferenceBasesLength(_) => f.write_str("invalid reference base length"),
            Self::PositionOverflow(start, len) => write!(
                f,
                "calculation of the end position overflowed: {start} + {len}",
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
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("ACGT".parse()?)
    ///     .set_info("END=8".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.end(), Ok(Position::from(8)));
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
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("ACGT".parse()?)
    ///     .build()?;
    ///
    /// assert_eq!(record.end(), Ok(Position::from(4)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> Result<Position, EndError> {
        use self::info::field::Value;
        use super::header::info::key;

        let end = if let Some(Some(value)) = self.info().get(&key::END_POSITION) {
            match value {
                Value::Integer(n) => usize::try_from(*n).map_err(EndError::InvalidPosition)?,
                _ => return Err(EndError::InvalidInfoEndPositionFieldValue),
            }
        } else {
            let start = usize::from(self.position());

            // `len` is guaranteed to be > 0.
            let len = self.reference_bases().len();

            start
                .checked_add(len - 1)
                .ok_or(EndError::PositionOverflow(start, len))?
        };

        Ok(Position::from(end))
    }
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{chrom}\t{pos}",
            chrom = self.chromosome(),
            pos = self.position()
        )?;

        if self.ids().is_empty() {
            write!(f, "\t{MISSING_FIELD}")?;
        } else {
            write!(f, "\t{id}", id = self.ids())?;
        }

        write!(f, "\t{ref}", r#ref = self.reference_bases())?;

        if self.alternate_bases().is_empty() {
            write!(f, "\t{MISSING_FIELD}")?;
        } else {
            write!(f, "\t{alt}", alt = self.alternate_bases())?;
        }

        if let Some(qual) = self.quality_score() {
            write!(f, "\t{qual}")?;
        } else {
            write!(f, "\t{MISSING_FIELD}")?;
        }

        if let Some(filter) = self.filters() {
            write!(f, "\t{filter}")?;
        } else {
            write!(f, "\t{MISSING_FIELD}")?;
        }

        if self.info().is_empty() {
            write!(f, "\t{MISSING_FIELD}")?;
        } else {
            write!(f, "\t{info}", info = self.info())?;
        }

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
    fn test_default() -> Result<(), Box<dyn std::error::Error>> {
        let actual = Record::default();

        let expected = Record::builder()
            .set_chromosome(".".parse()?)
            .set_position(Position::from(1))
            .set_reference_bases("N".parse()?)
            .build()?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_end() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::info::key;

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .set_info([(key::END_POSITION, None)].into_iter().collect())
            .build()?;

        assert_eq!(record.end(), Ok(Position::from(1)));

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .set_info(
                [(key::END_POSITION, Some(info::field::Value::Flag))]
                    .into_iter()
                    .collect(),
            )
            .build()?;

        assert_eq!(
            record.end(),
            Err(EndError::InvalidInfoEndPositionFieldValue)
        );

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::from(usize::MAX))
            .set_reference_bases("ACGT".parse()?)
            .build()?;

        assert_eq!(record.end(), Err(EndError::PositionOverflow(usize::MAX, 4)));

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .build()?;

        assert_eq!(record.to_string(), "sq0\t1\t.\tA\t.\t.\t.\t.");

        Ok(())
    }

    #[test]
    fn test_fmt_with_format() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .set_genotypes("GT:GQ\t0|0:13".parse()?)
            .build()?;

        assert_eq!(
            record.to_string(),
            "sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13"
        );

        Ok(())
    }
}
