use std::io;

use super::RecordBuf;
use crate::{record_buf::Attributes, Record};

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

        builder = builder.set_strand(record.strand()?);

        if let Some(frame) = record.frame().transpose()? {
            builder = builder.set_frame(frame);
        }

        let mut raw_attributes = Vec::new();

        for result in record.attributes()?.iter() {
            let (key, value) = result?;

            for v in value.iter() {
                raw_attributes.push((key.into(), v.into()));
            }
        }

        let attributes = Attributes::from(raw_attributes);
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
            .set_attributes(
                [(String::from("id"), String::from("0"))]
                    .into_iter()
                    .collect(),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
