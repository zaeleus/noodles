use std::io::{self, BufRead};

use super::read_line;
use crate::Record;

pub(super) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.buf.clear();

    let mut len = 0;

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.name_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.flags_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.reference_sequence_name_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.alignment_start_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.mapping_quality_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.cigar_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.mate_reference_sequence_name_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.mate_alignment_start_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.template_length_end = record.buf.len();

    len += read_required_field(reader, &mut record.buf)?;
    record.bounds.sequence_end = record.buf.len();

    let (n, is_eol) = read_last_required_field(reader, &mut record.buf)?;
    len += n;
    record.bounds.quality_scores_end = record.buf.len();

    if !is_eol {
        len += read_line(reader, &mut record.buf)?;
    }

    Ok(len)
}

fn read_required_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    let (len, is_eol) = read_field(reader, dst)?;

    if is_eol {
        Err(io::Error::new(io::ErrorKind::InvalidData, "unexpected EOL"))
    } else {
        Ok(len)
    }
}

fn read_last_required_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<(usize, bool)>
where
    R: BufRead,
{
    read_field(reader, dst)
}

fn read_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<(usize, bool)>
where
    R: BufRead,
{
    use memchr::memchr2;

    const DELIMITER: u8 = b'\t';
    const LINE_FEED: u8 = b'\n';

    let mut r#match = None;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if r#match.is_some() || src.is_empty() {
            break;
        }

        let n = match memchr2(DELIMITER, LINE_FEED, src) {
            Some(i) => {
                dst.extend(&src[..i]);
                r#match = Some(src[i]);
                i + 1
            }
            None => {
                dst.extend(src);
                src.len()
            }
        };

        len += n;

        reader.consume(n);
    }

    let is_eol = matches!(r#match, Some(LINE_FEED));

    Ok((len, is_eol))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"[..];

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

        let mut src = &b"\n"[..];
        assert!(matches!(
            read_record(&mut src, &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
