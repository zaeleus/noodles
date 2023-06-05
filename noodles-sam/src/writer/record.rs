mod cigar;
mod data;
mod position;
mod quality_scores;
mod sequence;

pub use self::{
    cigar::write_cigar, data::write_data, position::write_position,
    quality_scores::write_quality_scores, sequence::write_sequence,
};

use std::io::{self, Write};

use crate::{alignment::Record, Header};

const MISSING: u8 = b'*';

pub fn write_record<W>(writer: &mut W, header: &Header, record: &Record) -> io::Result<()>
where
    W: Write,
{
    use super::num;

    const DELIMITER: &[u8] = b"\t";
    const EQ: &[u8] = b"=";
    const MISSING: &[u8] = b"*";

    let qname = record
        .read_name()
        .map(|name| name.as_ref())
        .unwrap_or(MISSING);

    let reference_sequence = record.reference_sequence(header).transpose()?;

    let rname = reference_sequence
        .map(|(name, _)| name.as_bytes())
        .unwrap_or(MISSING);

    let mapq = record
        .mapping_quality()
        .map(u8::from)
        .unwrap_or(crate::record::mapping_quality::MISSING);

    let rnext = record
        .mate_reference_sequence(header)
        .transpose()?
        .map(|(mate_reference_sequence_name, _)| {
            if let Some((reference_sequence_name, _)) = reference_sequence {
                if mate_reference_sequence_name == reference_sequence_name {
                    return EQ;
                }
            }

            mate_reference_sequence_name.as_ref()
        })
        .unwrap_or(MISSING);

    writer.write_all(qname)?;

    writer.write_all(DELIMITER)?;
    num::write_u16(writer, u16::from(record.flags()))?;

    writer.write_all(DELIMITER)?;
    writer.write_all(rname)?;

    writer.write_all(DELIMITER)?;
    write_position(writer, record.alignment_start())?;

    writer.write_all(DELIMITER)?;
    num::write_u8(writer, mapq)?;

    writer.write_all(DELIMITER)?;
    write_cigar(writer, record.cigar())?;

    writer.write_all(DELIMITER)?;
    writer.write_all(rnext)?;

    writer.write_all(DELIMITER)?;
    write_position(writer, record.mate_alignment_start())?;

    writer.write_all(DELIMITER)?;
    num::write_i32(writer, record.template_length())?;

    writer.write_all(DELIMITER)?;
    write_sequence(writer, record.cigar().read_length(), record.sequence())?;

    writer.write_all(DELIMITER)?;
    write_quality_scores(writer, record.sequence().len(), record.quality_scores())?;

    write_data(writer, record.data())?;

    writeln!(writer)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record_with_data() -> io::Result<()> {
        use crate::record::data::field::{tag, Value};

        let mut buf = Vec::new();

        let header = Header::default();

        let data = [(tag::READ_GROUP, Value::String(String::from("rg0")))]
            .into_iter()
            .collect();
        let record = Record::builder().set_data(data).build();

        write_record(&mut buf, &header, &record)?;

        let expected = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
