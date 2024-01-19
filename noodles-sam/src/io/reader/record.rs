use std::io::{self, BufRead};

use super::read_line;
use crate::Record;

pub(super) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.buf.clear();

    let mut len = 0;

    len += read_field(reader, &mut record.buf)?;
    record.bounds.name_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.flags_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.reference_sequence_name_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.alignment_start_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.mapping_quality_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.cigar_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.mate_reference_sequence_name_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.mate_alignment_start_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.template_length_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.sequence_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.quality_scores_end = record.buf.len();

    len += read_line(reader, &mut record.buf)?;

    Ok(len)
}

fn read_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const DELIMITER: u8 = b'\t';

    let mut is_delimiter = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if is_delimiter || src.is_empty() {
            break;
        }

        let n = match src.iter().position(|&b| b == DELIMITER) {
            Some(i) => {
                dst.extend_from_slice(&src[..i]);
                is_delimiter = true;
                i + 1
            }
            None => {
                dst.extend_from_slice(src);
                src.len()
            }
        };

        len += n;

        reader.consume(n);
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*"[..];

        let mut record = Record::default();
        read_record(&mut src, &mut record)?;

        assert_eq!(record.buf, b"*4*0255**00**");

        assert_eq!(record.bounds.name_end, 1);
        assert_eq!(record.bounds.flags_end, 2);
        assert_eq!(record.bounds.reference_sequence_name_end, 3);
        assert_eq!(record.bounds.alignment_start_end, 4);
        assert_eq!(record.bounds.mapping_quality_end, 7);
        assert_eq!(record.bounds.cigar_end, 8);
        assert_eq!(record.bounds.mate_reference_sequence_name_end, 9);
        assert_eq!(record.bounds.mate_alignment_start_end, 10);
        assert_eq!(record.bounds.template_length_end, 11);
        assert_eq!(record.bounds.sequence_end, 12);
        assert_eq!(record.bounds.quality_scores_end, 13);

        Ok(())
    }
}
