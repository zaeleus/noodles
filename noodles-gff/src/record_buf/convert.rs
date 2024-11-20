use std::io;

use super::{attributes::field::Value, RecordBuf, MISSING_FIELD};
use crate::{record::attributes::field::Value as ValueRef, Record};

impl<'l> TryFrom<Record<'l>> for RecordBuf {
    type Error = io::Error;

    fn try_from(record: Record<'l>) -> Result<Self, Self::Error> {
        let mut builder = Self::builder();

        builder = builder
            .set_reference_sequence_name(record.reference_sequence_name().into())
            .set_source(record.source().into())
            .set_type(record.ty().into())
            .set_start(record.start()?)
            .set_end(record.end()?);

        if record.score() != MISSING_FIELD {
            let score = record
                .score()
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

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
                result.map(|(k, v)| {
                    let value = match v {
                        ValueRef::String(s) => Value::from(s),
                        ValueRef::Array(values) => {
                            Value::Array(values.iter().map(String::from).collect())
                        }
                    };

                    (k.into(), value)
                })
            })
            .collect::<io::Result<_>>()?;

        builder = builder.set_attributes(attributes);

        Ok(builder.build())
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::record_buf::Strand;

    #[test]
    fn test_try_from_record_for_record_buf() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::try_new("sq0\t.\texon\t8\t13\t.\t+\t.\tID=0;Name=n0")?;
        let actual = RecordBuf::try_from(record)?;

        let expected = RecordBuf::builder()
            .set_reference_sequence_name(String::from("sq0"))
            .set_source(String::from("."))
            .set_type(String::from("exon"))
            .set_start(Position::try_from(8)?)
            .set_end(Position::try_from(13)?)
            .set_strand(Strand::Forward)
            .set_attributes(
                [
                    (String::from("ID"), Value::from("0")),
                    (String::from("Name"), Value::from("n0")),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
