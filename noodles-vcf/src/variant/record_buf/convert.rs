use std::io;

use super::RecordBuf;
use crate::{Header, variant::Record};

impl RecordBuf {
    /// Converts a variant record to a buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::RecordBuf};
    ///
    /// let header = vcf::Header::default();
    /// let record = vcf::Record::default();
    ///
    /// let record_buf = RecordBuf::try_from_variant_record(&header, &record)?;
    ///
    /// assert_eq!(record_buf, RecordBuf::default());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn try_from_variant_record<R>(header: &Header, record: &R) -> io::Result<Self>
    where
        R: Record + ?Sized,
    {
        let mut dst = RecordBuf::default();
        dst.try_clone_from_variant_record(header, record)?;
        Ok(dst)
    }

    /// Clones the given record into this record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, variant::RecordBuf};
    ///
    /// let mut record_buf = RecordBuf::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .build();
    ///
    /// let header = vcf::Header::default();
    /// let record = vcf::Record::default();
    ///
    /// record_buf.try_clone_from_variant_record(&header, &record)?;
    ///
    /// assert_eq!(record_buf, RecordBuf::default());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn try_clone_from_variant_record<R>(
        &mut self,
        header: &Header,
        record: &R,
    ) -> io::Result<()>
    where
        R: Record + ?Sized,
    {
        let src_reference_sequence_name = record.reference_sequence_name(header)?;
        let dst_reference_sequence_name = self.reference_sequence_name_mut();
        dst_reference_sequence_name.clear();
        dst_reference_sequence_name.push_str(src_reference_sequence_name);

        *self.variant_start_mut() = record.variant_start().transpose()?;

        let dst_raw_ids = self.ids_mut().as_mut();
        dst_raw_ids.clear();
        dst_raw_ids.extend(record.ids().iter().map(String::from));

        let dst_reference_bases = self.reference_bases_mut();
        dst_reference_bases.clear();

        for result in record.reference_bases().iter() {
            let base = result?;
            dst_reference_bases.push(char::from(base));
        }

        let dst_alternate_bases = self.alternate_bases_mut().as_mut();
        dst_alternate_bases.clear();

        for result in record.alternate_bases().iter() {
            let base = result?;
            dst_alternate_bases.push(base.map(String::from));
        }

        *self.quality_score_mut() = record.quality_score().transpose()?;

        let dst_filters = self.filters_mut().as_mut();
        dst_filters.clear();

        for result in record.filters().iter(header) {
            let filter = result?;
            dst_filters.insert(filter.into());
        }

        let dst_info = self.info_mut().as_mut();
        dst_info.clear();

        for result in record.info().iter(header) {
            let (key, value) = result?;
            let v = value.map(|v| v.try_into()).transpose()?;
            dst_info.insert(key.into(), v);
        }

        let src_samples = record.samples()?;
        let dst_samples = self.samples_mut();

        let dst_samples_keys = dst_samples.keys_mut().as_mut();
        dst_samples_keys.clear();

        for result in src_samples.column_names(header) {
            let key = result?;
            dst_samples_keys.insert(key.into());
        }

        let dst_samples_values = &mut dst_samples.values;
        dst_samples_values.clear();

        for src_sample in src_samples.iter() {
            let dst_sample = src_sample
                .iter(header)
                .map(|result| result.and_then(|(_, value)| value.map(|v| v.try_into()).transpose()))
                .collect::<io::Result<_>>()?;

            dst_samples_values.push(dst_sample);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::variant::{
        record::{info, samples},
        record_buf::{
            Filters, Samples, info::field::Value as InfoValueBuf,
            samples::sample::Value as SampleValueBuf,
        },
    };

    #[test]
    fn test_try_clone_from_variant_record() -> io::Result<()> {
        let header = Header::builder().add_sample_name("s0").build();

        let keys = [String::from(samples::keys::key::GENOTYPE)]
            .into_iter()
            .collect();

        let samples = Samples::new(keys, vec![vec![Some(SampleValueBuf::from("0|0"))]]);

        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_variant_start(Position::MIN)
            .set_ids([String::from("nd0")].into_iter().collect())
            .set_reference_bases("A")
            .set_alternate_bases([Some(String::from("C"))].into_iter().collect())
            .set_filters(Filters::pass())
            .set_info(
                [(
                    String::from(info::field::key::SAMPLES_WITH_DATA_COUNT),
                    Some(InfoValueBuf::Integer(3)),
                )]
                .into_iter()
                .collect(),
            )
            .set_samples(samples)
            .build();

        let mut dst = RecordBuf::default();
        dst.try_clone_from_variant_record(&header, &record)?;

        assert_eq!(dst, record);

        Ok(())
    }
}
