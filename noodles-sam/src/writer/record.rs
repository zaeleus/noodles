mod cigar;
mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod sequence;
mod template_length;

pub use self::{
    cigar::write_cigar, data::write_data, position::write_position,
    quality_scores::write_quality_scores, sequence::write_sequence,
};

use std::io::{self, Write};

use noodles_core::Position;

use self::{
    flags::write_flags, mapping_quality::write_mapping_quality, name::write_name,
    template_length::write_template_length,
};
use crate::{
    alignment::{record_buf::Flags, Record},
    record::MappingQuality,
    Header,
};

const MISSING: u8 = b'*';

pub fn write_record<W, R>(writer: &mut W, header: &Header, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record + ?Sized,
{
    const DELIMITER: &[u8] = b"\t";
    const EQ: &[u8] = b"=";
    const MISSING: &[u8] = b"*";

    let reference_sequence = record.reference_sequence(header).transpose()?;
    let reference_sequence_name = reference_sequence.map(|(name, _)| name).unwrap_or(MISSING);

    let mate_reference_sequence_name = record
        .mate_reference_sequence(header)
        .transpose()?
        .map(|(mate_reference_sequence_name, _)| {
            if let Some((reference_sequence_name, _)) = reference_sequence {
                if mate_reference_sequence_name == reference_sequence_name {
                    return EQ;
                }
            }

            mate_reference_sequence_name
        })
        .unwrap_or(MISSING);

    write_name(writer, record.name())?;

    writer.write_all(DELIMITER)?;
    let flags = Flags::try_from(record.flags().as_ref())?;
    write_flags(writer, flags)?;

    writer.write_all(DELIMITER)?;
    writer.write_all(reference_sequence_name)?;

    writer.write_all(DELIMITER)?;
    let alignment_start = record
        .alignment_start()
        .map(|position| Position::try_from(position.as_ref()))
        .transpose()?;
    write_position(writer, alignment_start)?;

    writer.write_all(DELIMITER)?;

    let mapping_quality = record
        .mapping_quality()
        .map(|mapping_quality| MappingQuality::try_from(mapping_quality.as_ref()))
        .transpose()?;

    write_mapping_quality(writer, mapping_quality)?;

    let cigar = record.cigar();

    writer.write_all(DELIMITER)?;
    write_cigar(writer, &cigar)?;

    writer.write_all(DELIMITER)?;
    writer.write_all(mate_reference_sequence_name)?;

    writer.write_all(DELIMITER)?;
    let mate_alignment_start = record
        .mate_alignment_start()
        .map(|position| Position::try_from(position.as_ref()))
        .transpose()?;
    write_position(writer, mate_alignment_start)?;

    writer.write_all(DELIMITER)?;
    let template_length = i32::try_from(record.template_length().as_ref())?;
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
        use crate::alignment::{record::data::field::tag, record_buf::data::field::Value};

        let mut buf = Vec::new();

        let header = Header::default();

        let data = [(tag::READ_GROUP, Value::from("rg0"))]
            .into_iter()
            .collect();
        let record = RecordBuf::builder().set_data(data).build();

        write_record(&mut buf, &header, &record)?;

        let expected = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
