//! Variant record buffer.

pub mod alternate_bases;
pub mod builder;
pub mod filters;
pub mod ids;
pub mod info;
pub mod samples;

use noodles_core::Position;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, filters::Filters, ids::Ids, info::Info,
    samples::Samples,
};

use std::{error, fmt, num};

pub(crate) const MISSING_FIELD: &str = ".";

/// A VCF record.
///
/// A VCF record has 8 required fields: chromosome (`CHROM`), position (`POS`), IDs (`ID`),
/// reference bases (`REF`), alternate bases (`ALT`), quality score (`QUAL`), filters (`FILTER`),
/// and information (`INFO`).
///
/// Additionally, each record can have genotype information. This adds the extra `FORMAT` field and
/// a number of genotype fields.
#[derive(Clone, Debug, PartialEq)]
pub struct RecordBuf {
    chromosome: String,
    position: Option<Position>,
    ids: Ids,
    reference_bases: String,
    alternate_bases: AlternateBases,
    quality_score: Option<f32>,
    filters: Option<Filters>,
    info: Info,
    samples: Samples,
}

impl RecordBuf {
    /// Returns a builder to create a record from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let builder = vcf::variant::RecordBuf::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the chromosome of the record.
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
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    /// Returns a mutable reference to the chromosome.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let mut record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .build();
    ///
    /// *record.chromosome_mut() = String::from("sq1");
    ///
    /// assert_eq!(record.chromosome(), "sq1");
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
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_position(Position::MIN)
    ///     .build();
    ///
    /// assert_eq!(record.position(), Some(Position::MIN));
    /// ```
    pub fn position(&self) -> Option<Position> {
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
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// *record.position_mut() = Some(Position::MIN);
    /// assert_eq!(record.position(), Some(Position::MIN));
    /// ```
    pub fn position_mut(&mut self) -> &mut Option<Position> {
        &mut self.position
    }

    /// Returns a list of IDs of the record.
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
    pub fn ids(&self) -> &Ids {
        &self.ids
    }

    /// Returns a mutable reference to the IDs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Ids};
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    ///
    /// let ids: Ids = [String::from("nd0")].into_iter().collect();
    /// *record.ids_mut() = ids.clone();
    ///
    /// assert_eq!(record.ids(), &ids);
    /// ```
    pub fn ids_mut(&mut self) -> &mut Ids {
        &mut self.ids
    }

    /// Returns the reference bases of the record.
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
    pub fn reference_bases(&self) -> &str {
        &self.reference_bases
    }

    /// Returns a mutable reference to the reference bases of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let mut record = vcf::variant::RecordBuf::builder()
    ///     .set_reference_bases("A")
    ///     .build();
    ///
    /// *record.reference_bases_mut() = String::from("T");
    ///
    /// assert_eq!(record.reference_bases(), "T");
    /// ```
    pub fn reference_bases_mut(&mut self) -> &mut String {
        &mut self.reference_bases
    }

    /// Returns the alternate bases of the record.
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
    pub fn alternate_bases(&self) -> &AlternateBases {
        &self.alternate_bases
    }

    /// Returns a mutable reference to the alternate bases of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::AlternateBases};
    ///
    /// let mut record = vcf::variant::RecordBuf::builder()
    ///     .set_reference_bases("A")
    ///     .build();
    ///
    /// let alternate_bases = AlternateBases::from(vec![String::from("C")]);
    /// *record.alternate_bases_mut() = alternate_bases.clone();
    ///
    /// assert_eq!(record.alternate_bases(), &alternate_bases);
    /// ```
    pub fn alternate_bases_mut(&mut self) -> &mut AlternateBases {
        &mut self.alternate_bases
    }

    /// Returns the quality score of the record.
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
    pub fn quality_score(&self) -> Option<f32> {
        self.quality_score
    }

    /// Returns a mutable reference to the quality score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// *record.quality_score_mut() = Some(13.0);
    /// assert_eq!(record.quality_score(), Some(13.0));
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
    /// use noodles_vcf::{self as vcf, variant::record_buf::Filters};
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_filters(Filters::Pass)
    ///     .build();
    ///
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
    /// ```
    pub fn filters(&self) -> Option<&Filters> {
        self.filters.as_ref()
    }

    /// Returns a mutable reference to the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Filters};
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// *record.filters_mut() = Some(Filters::Pass);
    /// assert_eq!(record.filters(), Some(&Filters::Pass));
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
    ///     variant::record_buf::{
    ///         info::field::{key, Value},
    ///         Info,
    ///     },
    /// };
    ///
    /// let info: Info = [
    ///     (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::from(3))),
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
    ///     variant::record_buf::{
    ///         info::field::{key, Value},
    ///         Info,
    ///     },
    /// };
    ///
    /// let info: Info = [
    ///     (String::from(key::SAMPLES_WITH_DATA_COUNT), Some(Value::from(3))),
    ///     (String::from(key::ALLELE_FREQUENCIES), Some(Value::from(vec![Some(0.5)]))),
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// let mut record = vcf::variant::RecordBuf::builder()
    ///     .set_info(info)
    ///     .build();
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
    /// let samples = Samples::new(
    ///     keys.clone(),
    ///     vec![vec![Some(Value::from("0|0")), Some(Value::from(13))]],
    /// );
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_samples(samples)
    ///     .build();
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
    ///     variant::record_buf::{
    ///         samples::{keys::key, sample::Value, Keys},
    ///         Samples
    ///     },
    /// };
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
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

impl Default for RecordBuf {
    fn default() -> Self {
        Self {
            chromosome: String::from("."),
            position: Some(Position::MIN),
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
    PositionOverflow(Position, usize),
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

impl RecordBuf {
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
    /// use noodles_core::Position;
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     variant::record_buf::info::field::{key, Value},
    /// };
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::MIN)
    ///     .set_reference_bases("ACGT")
    ///     .set_info(
    ///         [(String::from(key::END_POSITION), Some(Value::from(8)))]
    ///             .into_iter()
    ///             .collect(),
    ///     )
    ///     .build();
    ///
    /// assert_eq!(record.end(), Ok(Position::try_from(8)?));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    ///
    /// ## Calculated using the start position and reference bases length
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::MIN)
    ///     .set_reference_bases("ACGT")
    ///     .build();
    ///
    /// assert_eq!(record.end(), Ok(Position::try_from(4)?));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn end(&self) -> Result<Position, EndError> {
        use self::info::field::{key, Value};

        if let Some(Some(value)) = self.info().get(key::END_POSITION) {
            match value {
                Value::Integer(n) => usize::try_from(*n)
                    .and_then(Position::try_from)
                    .map_err(EndError::InvalidPosition),
                _ => return Err(EndError::InvalidInfoEndPositionFieldValue),
            }
        } else {
            let start = self.position().unwrap_or(Position::MIN);
            let reference_bases = self.reference_bases();

            let len = if reference_bases.is_empty() {
                return Err(EndError::InvalidReferenceBasesLength { actual: 0 });
            } else {
                reference_bases.len()
            };

            start
                .checked_add(len - 1)
                .ok_or(EndError::PositionOverflow(start, len))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() -> Result<(), Box<dyn std::error::Error>> {
        let actual = RecordBuf::default();

        let expected = RecordBuf::builder()
            .set_chromosome(".")
            .set_position(Position::MIN)
            .set_reference_bases("N")
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_end() -> Result<(), Box<dyn std::error::Error>> {
        use super::info::field::key;

        let record = RecordBuf::builder()
            .set_chromosome("sq0")
            .set_position(Position::MIN)
            .set_reference_bases("A")
            .set_info(
                [(String::from(key::END_POSITION), None)]
                    .into_iter()
                    .collect(),
            )
            .build();

        assert_eq!(record.end(), Ok(Position::MIN));

        let record = RecordBuf::builder()
            .set_chromosome("sq0")
            .set_position(Position::MIN)
            .set_reference_bases("A")
            .set_info(
                [(
                    String::from(key::END_POSITION),
                    Some(info::field::Value::Flag),
                )]
                .into_iter()
                .collect(),
            )
            .build();

        assert_eq!(
            record.end(),
            Err(EndError::InvalidInfoEndPositionFieldValue)
        );

        let record = RecordBuf::builder()
            .set_chromosome("sq0")
            .set_position(Position::MAX)
            .set_reference_bases("ACGT")
            .build();

        assert_eq!(
            record.end(),
            Err(EndError::PositionOverflow(Position::MAX, 4))
        );

        Ok(())
    }
}
