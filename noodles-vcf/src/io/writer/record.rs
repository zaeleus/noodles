mod alternate_bases;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_sequence_name;
mod samples;

use std::io::{self, Write};

use self::{
    alternate_bases::write_alternate_bases, filters::write_filters, ids::write_ids,
    info::write_info, position::write_position, quality_score::write_quality_score,
    reference_sequence_name::write_reference_sequence_name, samples::write_samples,
};
use crate::{variant::RecordBuf, Header};

const MISSING: &[u8] = b".";

pub(super) fn write_record<W>(writer: &mut W, header: &Header, record: &RecordBuf) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b"\t";

    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    writer.write_all(DELIMITER)?;
    write_position(writer, record.position())?;

    writer.write_all(DELIMITER)?;
    write_ids(writer, record.ids())?;

    writer.write_all(DELIMITER)?;
    write!(writer, "{}", record.reference_bases())?;

    writer.write_all(DELIMITER)?;
    write_alternate_bases(writer, record.alternate_bases())?;

    writer.write_all(DELIMITER)?;
    write_quality_score(writer, record.quality_score())?;

    writer.write_all(DELIMITER)?;
    write_filters(writer, record.filters())?;

    writer.write_all(DELIMITER)?;
    write_info(writer, header, record.info())?;

    if !record.samples().is_empty() {
        writer.write_all(DELIMITER)?;
        write_samples(writer, record.samples())?;
    }

    writer.write_all(b"\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_position(Position::MIN)
            .set_reference_bases("A")
            .build();

        let header = Header::default();
        let mut buf = Vec::new();
        write_record(&mut buf, &header, &record)?;
        assert_eq!(buf, b"sq0\t1\t.\tA\t.\t.\t.\t.\n");

        Ok(())
    }
}
