use std::{io, str};

use super::{attributes::field::Value, RecordBuf};
use crate::feature::{record::attributes::field::Value as ValueRef, Record};

impl RecordBuf {
    /// Converts a feature record to a record buffer.
    pub fn try_from_feature_record<R>(record: &R) -> io::Result<Self>
    where
        R: Record,
    {
        let mut builder = Self::builder();

        builder = builder
            .set_reference_sequence_name(str::from_utf8(record.reference_sequence_name()).unwrap()) // FIXME
            .set_source(str::from_utf8(record.source()).unwrap()) // FIXME
            .set_type(str::from_utf8(record.ty()).unwrap()) //FIXME
            .set_start(record.feature_start()?)
            .set_end(record.feature_end()?);

        if let Some(score) = record.score().transpose()? {
            builder = builder.set_score(score);
        }

        builder = builder.set_strand(record.strand()?);

        if let Some(phase) = record.phase().transpose()? {
            builder = builder.set_phase(phase);
        }

        let attributes = record
            .attributes()
            .iter()
            .map(|result| {
                result.and_then(|(k, v)| {
                    let value = match v {
                        ValueRef::String(s) => Value::String(s.into_owned()),
                        ValueRef::Array(values) => Value::Array(
                            values
                                .iter()
                                .map(|result| result.map(|s| s.into_owned()))
                                .collect::<io::Result<_>>()?,
                        ),
                    };

                    Ok((k.into_owned(), value))
                })
            })
            .collect::<io::Result<_>>()?;

        builder = builder.set_attributes(attributes);

        Ok(builder.build())
    }
}
