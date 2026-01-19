mod alternate_bases;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_bases;
mod reference_sequence_name;
mod samples;

use std::io::{self, Write};

use self::{
    alternate_bases::write_alternate_bases, filters::write_filters, ids::write_ids,
    info::write_info, position::write_position, quality_score::write_quality_score,
    reference_bases::write_reference_bases, reference_sequence_name::write_reference_sequence_name,
    samples::write_samples,
};
use crate::{Header, variant::Record};

const MISSING: &[u8] = b".";

pub(super) fn write_record<W, R>(writer: &mut W, header: &Header, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record + ?Sized,
{
    const DELIMITER: &[u8] = b"\t";

    let reference_sequence_name = record.reference_sequence_name(header)?;
    write_reference_sequence_name(writer, reference_sequence_name)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_all(DELIMITER)?;
    let position = record.variant_start().transpose()?;
    write_position(writer, position)?;

    writer.write_all(DELIMITER)?;
    write_ids(writer, record.ids()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_all(DELIMITER)?;
    write_reference_bases(writer, record.reference_bases())?;

    writer.write_all(DELIMITER)?;
    write_alternate_bases(writer, record.alternate_bases())?;

    writer.write_all(DELIMITER)?;
    let quality_score = record.quality_score().transpose()?;
    write_quality_score(writer, quality_score)?;

    writer.write_all(DELIMITER)?;
    write_filters(writer, header, record.filters())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_all(DELIMITER)?;
    write_info(writer, header, record.info())?;

    let samples = record.samples()?;

    if !samples.is_empty() {
        writer.write_all(DELIMITER)?;
        write_samples(writer, header, samples)?;
    }

    writer.write_all(b"\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::variant::RecordBuf;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_variant_start(Position::MIN)
            .set_reference_bases("A")
            .build();

        let header = Header::default();
        let mut buf = Vec::new();
        write_record(&mut buf, &header, &record)?;
        assert_eq!(buf, b"sq0\t1\t.\tA\t.\t.\t.\t.\n");

        Ok(())
    }
}
