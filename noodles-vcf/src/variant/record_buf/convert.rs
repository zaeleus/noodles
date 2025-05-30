use std::io;

use super::RecordBuf;
use crate::{Header, variant::Record};

impl RecordBuf {
    /// Converts a variant record to a buffer.
    pub fn try_from_variant_record<R>(header: &Header, record: &R) -> io::Result<Self>
    where
        R: Record + ?Sized,
    {
        use super::Samples;

        let mut record_buf = RecordBuf::default();

        *record_buf.reference_sequence_name_mut() = record.reference_sequence_name(header)?.into();
        *record_buf.variant_start_mut() = record.variant_start().transpose()?;
        *record_buf.ids_mut() = record.ids().iter().map(String::from).collect();

        let raw_reference_bases: Vec<_> =
            record.reference_bases().iter().collect::<io::Result<_>>()?;
        *record_buf.reference_bases_mut() = String::from_utf8(raw_reference_bases)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let raw_alternate_bases: Vec<_> = record
            .alternate_bases()
            .iter()
            .map(|result| result.map(String::from))
            .collect::<io::Result<_>>()?;
        *record_buf.alternate_bases_mut() = raw_alternate_bases.into();

        *record_buf.quality_score_mut() = record.quality_score().transpose()?;

        *record_buf.filters_mut() = record
            .filters()
            .iter(header)
            .map(|result| result.map(String::from))
            .collect::<io::Result<_>>()?;

        *record_buf.info_mut() = record
            .info()
            .iter(header)
            .map(|result| {
                result.and_then(|(key, value)| {
                    let v = value.map(|v| v.try_into()).transpose()?;
                    Ok((String::from(key), v))
                })
            })
            .collect::<io::Result<_>>()?;

        let samples = record.samples()?;

        let keys = samples
            .column_names(header)
            .map(|result| result.map(String::from))
            .collect::<io::Result<_>>()?;

        let values = samples
            .iter()
            .map(|sample| {
                sample
                    .iter(header)
                    .map(|result| {
                        result.and_then(|(_, value)| value.map(|v| v.try_into()).transpose())
                    })
                    .collect()
            })
            .collect::<io::Result<_>>()?;

        *record_buf.samples_mut() = Samples::new(keys, values);

        Ok(record_buf)
    }
}
