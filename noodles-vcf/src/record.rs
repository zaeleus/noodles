//! VCF record and fields.

pub mod alternate_bases;
pub mod builder;
pub mod filters;
pub mod ids;
pub mod info;
pub mod position;
pub mod samples;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, filters::Filters, ids::Ids, info::Info,
    position::Position, samples::Samples,
};

use std::{error, fmt, num};

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
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(record.chromosome(), "sq0");
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    /// Returns a mutable reference to the chromosome.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// *record.chromosome_mut() = String::from("sq1");
    ///
    /// assert_eq!(record.chromosome(), "sq1");
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn chromosome_mut(&mut self) -> &mut String {
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(8))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(usize::from(record.position()), 8);
    /// # Ok::<_, vcf::record::builder::BuildError>(())
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(8))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// *record.position_mut() = Position::from(13);
    ///
    /// assert_eq!(usize::from(record.position()), 13);
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn position_mut(&mut self) -> &mut Position {
        &mut self.position
    }

    /// Returns a list of IDs of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{Ids, Position}};
    ///
    /// let ids: Ids = [String::from("nd0")].into_iter().collect();
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_ids(ids.clone())
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(record.ids(), &ids);
    /// # Ok::<(), vcf::record::builder::BuildError>(())
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// let ids: Ids = [String::from("nd0")].into_iter().collect();
    /// *record.ids_mut() = ids.clone();
    ///
    /// assert_eq!(record.ids(), &ids);
    /// # Ok::<(), vcf::record::builder::BuildError>(())
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
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// assert_eq!(record.reference_bases(), "A");
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn reference_bases(&self) -> &str {
        &self.reference_bases
    }

    /// Returns a mutable reference to the reference bases of the record.
    ///
    /// This is a required field and guaranteed to be nonempty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// *record.reference_bases_mut() = String::from("T");
    ///
    /// assert_eq!(record.reference_bases(), "T");
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn reference_bases_mut(&mut self) -> &mut String {
        &mut self.reference_bases
    }

    /// Returns the alternate bases of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{AlternateBases, Position}};
    ///
    /// let alternate_bases = AlternateBases::from(vec![String::from("C")]);
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_alternate_bases(alternate_bases.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.alternate_bases(), &alternate_bases);
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn alternate_bases(&self) -> &AlternateBases {
        &self.alternate_bases
    }

    /// Returns a mutable reference to the alternate bases of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::{AlternateBases, Position}};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// let alternate_bases = AlternateBases::from(vec![String::from("C")]);
    /// *record.alternate_bases_mut() = alternate_bases.clone();
    ///
    /// assert_eq!(record.alternate_bases(), &alternate_bases);
    /// # Ok::<_, vcf::record::builder::BuildError>(())
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
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_quality_score(13.0)
    ///     .build()?;
    ///
    /// assert_eq!(record.quality_score(), Some(13.0));
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn quality_score(&self) -> Option<f32> {
        self.quality_score
    }

    /// Returns a mutable reference to the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// *record.quality_score_mut() = Some(13.0);
    ///
    /// assert_eq!(record.quality_score(), Some(13.0));
    /// # Ok::<_, vcf::record::builder::BuildError>(())
    /// ```
    pub fn quality_score_mut(&mut self) -> &mut Option<f32> {
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_filters(Filters::Pass)
    ///     .build()?;
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// # Ok::<_, vcf::record::builder::BuildError>(())
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
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// *record.filters_mut() = Some(Filters::Pass);
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// # Ok::<_, vcf::record::builder::BuildError>(())
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
    ///     record::{info::field::{key, Value}, Info, Position},
    /// };
    ///
    /// let info: Info = [
    ///     (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::from(3))),
    ///     (String::from(key::ALLELE_FREQUENCIES), Some(Value::from(vec![Some(0.5)]))),
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_info(info.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.info(), &info);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
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
    ///     record::{info::field::{key, Value}, Info, Position},
    /// };
    ///
    /// let info: Info = [
    ///     (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::from(3))),
    ///     (String::from(key::ALLELE_FREQUENCIES), Some(Value::from(vec![Some(0.5)]))),
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_info(info)
    ///     .build()?;
    ///
    /// record.info_mut().insert(String::from(key::TOTAL_DEPTH), Some(Value::Integer(13)));
    ///
    /// let expected = [
    ///     (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::Integer(3))),
    ///     (String::from(key::ALLELE_FREQUENCIES), Some(Value::from(vec![Some(0.5)]))),
    ///     (String::from(key::TOTAL_DEPTH), Some(Value::Integer(13))),
    /// ]
    /// .into_iter()
    /// .collect();
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
    ///     record::{
    ///         samples::{keys::key, sample::Value, Keys},
    ///         Position, Samples,
    ///     },
    /// };
    ///
    /// let keys = Keys::try_from(vec![
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ])?;
    /// let samples = Samples::new(
    ///     keys.clone(),
    ///     vec![vec![Some(Value::from("0|0")), Some(Value::from(13))]],
    /// );
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_samples(samples)
    ///     .build()?;
    ///
    /// assert_eq!(record.format(), &keys);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn format(&self) -> &samples::Keys {
        self.samples.keys()
    }

    /// Returns the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{
    ///         samples::{keys::key, sample::Value, Keys},
    ///         Position, Samples,
    ///     },
    /// };
    ///
    /// let keys = Keys::try_from(vec![
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ])?;
    /// let samples = Samples::new(
    ///     keys,
    ///     vec![vec![Some(Value::from("0|0")), Some(Value::from(13))]],
    /// );
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .set_samples(samples.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.samples(), &samples);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// Returns a mutable reference to the genotypes of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{
    ///         samples::{keys::key, sample::Value, Keys},
    ///         Position, Samples
    ///     },
    /// };
    ///
    /// let mut record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A")
    ///     .build()?;
    ///
    /// let keys = Keys::try_from(vec![
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ])?;
    /// let samples = Samples::new(
    ///     keys,
    ///     vec![vec![Some(Value::from("0|0")), Some(Value::from(13))]],
    /// );
    ///
    /// *record.samples_mut() = samples.clone();
    ///
    /// assert_eq!(record.samples(), &samples);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn samples_mut(&mut self) -> &mut Samples {
        &mut self.samples
    }
}

impl Default for Record {
    fn default() -> Self {
        Self {
            chromosome: String::from("."),
            position: Position::from(1),
            ids: Ids::default(),
            reference_bases: String::from("N"),
            alternate_bases: AlternateBases::default(),
            quality_score: None,
            filters: None,
            info: Info::default(),
            samples: Samples::default(),
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
    /// The reference bases length is invalid.
    InvalidReferenceBasesLength {
        /// The actual length.
        actual: usize,
    },
    /// The calculation of the end position overflowed.
    PositionOverflow(usize, usize),
}

impl error::Error for EndError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidPosition(e) => Some(e),
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
            Self::InvalidReferenceBasesLength { actual } => write!(
                f,
                "invalid reference base length: expected > 0, got {actual}"
            ),
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
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     record::{info::field::{key, Value}, Position},
    /// };
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("ACGT")
    ///     .set_info(
    ///         [(String::from(key::END_POSITION), Some(Value::from(8)))]
    ///             .into_iter()
    ///             .collect(),
    ///     )
    ///     .build()?;
    ///
    /// assert_eq!(record.end(), Ok(Position::from(8)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// ## Calculated using the start position and reference bases length
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(record.end(), Ok(Position::from(4)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> Result<Position, EndError> {
        use self::info::field::{key, Value};

        let end = if let Some(Some(value)) = self.info().get(key::END_POSITION) {
            match value {
                Value::Integer(n) => usize::try_from(*n).map_err(EndError::InvalidPosition)?,
                _ => return Err(EndError::InvalidInfoEndPositionFieldValue),
            }
        } else {
            let start = usize::from(self.position());
            let reference_bases = self.reference_bases();

            let len = if reference_bases.is_empty() {
                return Err(EndError::InvalidReferenceBasesLength { actual: 0 });
            } else {
                reference_bases.len()
            };

            start
                .checked_add(len - 1)
                .ok_or(EndError::PositionOverflow(start, len))?
        };

        Ok(Position::from(end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() -> Result<(), Box<dyn std::error::Error>> {
        let actual = Record::default();

        let expected = Record::builder()
            .set_chromosome(".")
            .set_position(Position::from(1))
            .set_reference_bases("N")
            .build()?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_end() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::info::field::key;

        let record = Record::builder()
            .set_chromosome("sq0")
            .set_position(Position::from(1))
            .set_reference_bases("A")
            .set_info(
                [(String::from(key::END_POSITION), None)]
                    .into_iter()
                    .collect(),
            )
            .build()?;

        assert_eq!(record.end(), Ok(Position::from(1)));

        let record = Record::builder()
            .set_chromosome("sq0")
            .set_position(Position::from(1))
            .set_reference_bases("A")
            .set_info(
                [(
                    String::from(key::END_POSITION),
                    Some(info::field::Value::Flag),
                )]
                .into_iter()
                .collect(),
            )
            .build()?;

        assert_eq!(
            record.end(),
            Err(EndError::InvalidInfoEndPositionFieldValue)
        );

        let record = Record::builder()
            .set_chromosome("sq0")
            .set_position(Position::from(usize::MAX))
            .set_reference_bases("ACGT")
            .build()?;

        assert_eq!(record.end(), Err(EndError::PositionOverflow(usize::MAX, 4)));

        Ok(())
    }
}
