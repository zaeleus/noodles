mod cigar;
mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_name;
mod sequence;
mod template_length;

pub use self::{
    cigar::write_cigar, data::write_data, position::write_position,
    quality_scores::write_quality_scores, sequence::write_sequence,
};

use std::io::{self, Write};

use self::{
    flags::write_flags,
    mapping_quality::write_mapping_quality,
    name::write_name,
    reference_sequence_name::{write_mate_reference_sequence_name, write_reference_sequence_name},
    template_length::write_template_length,
};
use crate::{alignment::Record, Header};

const MISSING: u8 = b'*';

pub fn write_record<W, R>(writer: &mut W, header: &Header, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record + ?Sized,
{
    const DELIMITER: &[u8] = b"\t";

    write_name(writer, record.name())?;

    writer.write_all(DELIMITER)?;
    let flags = record.flags()?;
    write_flags(writer, flags)?;

    writer.write_all(DELIMITER)?;

    let reference_sequence_name = record
        .reference_sequence(header)
        .transpose()?
        .map(|(name, _)| name.as_ref());

    write_reference_sequence_name(writer, reference_sequence_name)?;

    writer.write_all(DELIMITER)?;
    let alignment_start = record.alignment_start().transpose()?;
    write_position(writer, alignment_start)?;

    writer.write_all(DELIMITER)?;
    let mapping_quality = record.mapping_quality().transpose()?;
    write_mapping_quality(writer, mapping_quality)?;

    let cigar = record.cigar();

    writer.write_all(DELIMITER)?;
    write_cigar(writer, &cigar)?;

    writer.write_all(DELIMITER)?;

    let mate_reference_sequence_name = record
        .mate_reference_sequence(header)
        .transpose()?
        .map(|(name, _)| name.as_ref());

    write_mate_reference_sequence_name(
        writer,
        reference_sequence_name,
        mate_reference_sequence_name,
    )?;

    writer.write_all(DELIMITER)?;
    let mate_alignment_start = record.mate_alignment_start().transpose()?;
    write_position(writer, mate_alignment_start)?;

    writer.write_all(DELIMITER)?;
    let template_length = record.template_length()?;
    write_template_length(writer, template_length)?;

    let sequence = record.sequence();
    let base_count = sequence.len();

    writer.write_all(DELIMITER)?;
    let read_length = cigar.read_length()?;
    write_sequence(writer, read_length, sequence)?;

    writer.write_all(DELIMITER)?;
    write_quality_scores(writer, base_count, record.quality_scores())?;

    write_data(writer, record.data())?;

    writeln!(writer)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::RecordBuf;

    #[test]
    fn test_write_record_with_data() -> io::Result<()> {
        use crate::alignment::{record::data::field::Tag, record_buf::data::field::Value};

        let mut buf = Vec::new();

        let header = Header::default();

        let data = [(Tag::READ_GROUP, Value::from("rg0"))]
            .into_iter()
            .collect();
        let record = RecordBuf::builder().set_data(data).build();

        write_record(&mut buf, &header, &record)?;

        let expected = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
