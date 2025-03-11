use std::io;

use super::{Attributes, RecordBuf};
use crate::Record;

impl<'l> TryFrom<Record<'l>> for RecordBuf {
    type Error = io::Error;

    fn try_from(record: Record<'l>) -> Result<Self, Self::Error> {
        let mut builder = Self::builder();

        builder = builder
            .set_reference_sequence_name(record.reference_sequence_name())
            .set_source(record.source())
            .set_type(record.ty())
            .set_start(record.start()?)
            .set_end(record.end()?);

        if let Some(score) = record.score().transpose()? {
            builder = builder.set_score(score);
        }

        if let Some(strand) = record.strand().transpose()? {
            builder = builder.set_strand(strand);
        }

        if let Some(frame) = record.frame().transpose()? {
            builder = builder.set_frame(frame);
        }

        let fields: Vec<_> = record
            .attributes()
            .iter()
            .map(|result| result.map(|(k, v)| (k.into(), v.into())))
            .collect::<io::Result<_>>()?;

        let attributes = Attributes::from(fields);
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
    fn test_try_from_record_for_record_buf() -> io::Result<()> {
        const START: Position = Position::new(8).unwrap();
        const END: Position = Position::new(13).unwrap();

        let record = Record::try_new("sq0\tNDLS\texon\t8\t13\t.\t+\t.\tid 0;")?;
        let actual = RecordBuf::try_from(record)?;

        let expected = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_source("NDLS")
            .set_type("exon")
            .set_start(START)
            .set_end(END)
            .set_strand(Strand::Forward)
            .set_attributes(Attributes::from(vec![(
                String::from("id"),
                String::from("0"),
            )]))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
