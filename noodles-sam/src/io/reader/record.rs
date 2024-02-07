use std::io::{self, BufRead};

use super::read_line;
use crate::Record;

pub(crate) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    let fields = record.fields_mut();

    let buf = &mut fields.buf;
    buf.clear();

    let bounds = &mut fields.bounds;

    let mut len = 0;

    len += read_required_field(reader, buf)?;
    bounds.name_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.flags_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.reference_sequence_name_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.alignment_start_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.mapping_quality_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.cigar_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.mate_reference_sequence_name_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.mate_alignment_start_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.template_length_end = buf.len();

    len += read_required_field(reader, buf)?;
    bounds.sequence_end = buf.len();

    let (n, is_eol) = read_last_required_field(reader, buf)?;
    len += n;
    bounds.quality_scores_end = buf.len();

    if !is_eol {
        len += read_line(reader, buf)?;
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
    const CARRIAGE_RETURN: u8 = b'\r';

    let mut r#match = None;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if r#match.is_some() || src.is_empty() {
            break;
        }

        let (mut buf, n) = match memchr2(DELIMITER, LINE_FEED, src) {
            Some(i) => {
                r#match = Some(src[i]);
                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        if let [head @ .., CARRIAGE_RETURN] = buf {
            buf = head;
        }

        dst.extend(buf);
        len += n;
        reader.consume(n);
    }

    let is_eol = matches!(r#match, Some(LINE_FEED));

    Ok((len, is_eol))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::fields::Bounds;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"[..];
        let mut record = Record::default();
        read_record(&mut src, &mut record)?;
        assert_eq!(record.fields().buf, b"*4*0255**00**");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\r\n"[..];
        let mut record = Record::default();
        read_record(&mut src, &mut record)?;
        assert_eq!(record.fields().buf, b"*4*0255**00**");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"\n"[..];
        assert!(matches!(
            read_record(&mut src, &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
