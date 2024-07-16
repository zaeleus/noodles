//! Variant record buffer.

mod alternate_bases;
pub mod builder;
mod convert;
mod filters;
pub mod ids;
pub mod info;
pub mod samples;

use std::io;

use noodles_core::Position;

pub use self::{
    alternate_bases::AlternateBases, builder::Builder, filters::Filters, ids::Ids, info::Info,
    samples::Samples,
};
use crate::Header;

/// A variant record buffer.
#[derive(Clone, Debug, PartialEq)]
pub struct RecordBuf {
    reference_sequence_name: String,
    variant_start: Option<Position>,
    ids: Ids,
    reference_bases: String,
    alternate_bases: AlternateBases,
    quality_score: Option<f32>,
    filters: Filters,
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

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// ```
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    /// Returns a mutable reference to the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let mut record = vcf::variant::RecordBuf::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .build();
    ///
    /// *record.reference_sequence_name_mut() = String::from("sq1");
    ///
    /// assert_eq!(record.reference_sequence_name(), "sq1");
    /// ```
    pub fn reference_sequence_name_mut(&mut self) -> &mut String {
        &mut self.reference_sequence_name
    }

    /// Returns the variant start position.
    ///
    /// This position is 1-based, inclusive. If the record represents the start of a telomeric
    /// breakend, this returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_variant_start(Position::MIN)
    ///     .build();
    ///
    /// assert_eq!(record.variant_start(), Some(Position::MIN));
    /// ```
    pub fn variant_start(&self) -> Option<Position> {
        self.variant_start
    }

    /// Returns a mutable reference to the variant start position.
    ///
    /// This position is 1-based, inclusive. If the record represents the start of a telomeric
    /// breakend, this returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// *record.variant_start_mut() = Some(Position::MIN);
    /// assert_eq!(record.variant_start(), Some(Position::MIN));
    /// ```
    pub fn variant_start_mut(&mut self) -> &mut Option<Position> {
        &mut self.variant_start
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
    ///     .set_filters(Filters::pass())
    ///     .build();
    ///
    /// assert!(record.filters().is_pass());
    /// ```
    pub fn filters(&self) -> &Filters {
        &self.filters
    }

    /// Returns a mutable reference to the filters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::record_buf::Filters};
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// *record.filters_mut() = Filters::pass();
    /// assert!(record.filters().is_pass());
    /// ```
    pub fn filters_mut(&mut self) -> &mut Filters {
        &mut self.filters
    }

    /// Returns the addition information of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     variant::{
    ///         record::info::field::key,
    ///         record_buf::{info::field::Value, Info},
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
    ///     variant::{
    ///         record::info::field::key,
    ///         record_buf::{info::field::Value, Info},
    ///     }
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
    ///     variant::{
    ///         record::samples::keys::key,
    ///         record_buf::{samples::{sample::Value, Keys}, Samples},
    ///     },
    /// };
    ///
    /// let keys: Keys = [
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ].into_iter().collect();
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
    ///     variant::{
    ///         record::samples::keys::key,
    ///         record_buf::{samples::sample::Value, Samples},
    ///     },
    /// };
    ///
    /// let keys = [
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ].into_iter().collect();
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
    ///     variant::{
    ///         record::samples::keys::key,
    ///         record_buf::{samples::sample::Value, Samples},
    ///     },
    /// };
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    ///
    /// let keys = [
    ///     String::from(key::GENOTYPE),
    ///     String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
    /// ].into_iter().collect();
    /// let samples = Samples::new(
    ///     keys,
    ///     vec![vec![Some(Value::from("0|0")), Some(Value::from(13))]],
    /// );
    ///
    /// *record.samples_mut() = samples.clone();
    ///
    /// assert_eq!(record.samples(), &samples);
    /// ```
    pub fn samples_mut(&mut self) -> &mut Samples {
        &mut self.samples
    }
}

impl Default for RecordBuf {
    fn default() -> Self {
        Self {
            reference_sequence_name: String::from("."),
            variant_start: Some(Position::MIN),
            ids: Ids::default(),
            reference_bases: String::from("N"),
            alternate_bases: AlternateBases::default(),
            quality_score: None,
            filters: Filters::default(),
            info: Info::default(),
            samples: Samples::default(),
        }
    }
}

impl super::Record for RecordBuf {
    fn reference_sequence_name<'a, 'h: 'a>(&'a self, _: &'h Header) -> io::Result<&'a str> {
        Ok(self.reference_sequence_name())
    }

    fn variant_start(&self) -> Option<std::io::Result<Position>> {
        self.variant_start().map(Ok)
    }

    fn ids(&self) -> Box<dyn super::record::Ids + '_> {
        Box::new(self.ids())
    }

    fn reference_bases(&self) -> Box<dyn super::record::ReferenceBases + '_> {
        Box::new(self.reference_bases())
    }

    fn alternate_bases(&self) -> Box<dyn super::record::AlternateBases + '_> {
        Box::new(self.alternate_bases())
    }

    fn quality_score(&self) -> Option<std::io::Result<f32>> {
        self.quality_score().map(Ok)
    }

    fn filters(&self) -> Box<dyn super::record::Filters + '_> {
        Box::new(self.filters())
    }

    fn info(&self) -> Box<dyn super::record::Info + '_> {
        Box::new(self.info())
    }

    fn samples(&self) -> io::Result<Box<dyn super::record::Samples + '_>> {
        Ok(Box::new(self.samples()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() -> Result<(), Box<dyn std::error::Error>> {
        let actual = RecordBuf::default();

        let expected = RecordBuf::builder()
            .set_reference_sequence_name(".")
            .set_variant_start(Position::MIN)
            .set_reference_bases("N")
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
