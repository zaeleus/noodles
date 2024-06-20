use std::io::{self, BufRead};

use crate::Record;

pub(super) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    let fields = &mut record.0;

    let dst = &mut fields.buf;
    dst.clear();

    let bounds = &mut fields.bounds;
    bounds.other_fields_ends.clear();

    let mut len = 0;

    len += read_required_field(reader, dst)?;
    bounds.reference_sequence_name_end = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.feature_start_end = dst.len();

    let (n, is_eol) = read_field(reader, dst)?;
    len += n;
    bounds.feature_end_end = dst.len();

    if is_eol {
        return Ok(len);
    }

    loop {
        let (n, is_eol) = read_field(reader, dst)?;

        if n == 0 {
            break;
        }

        len += n;
        bounds.other_fields_ends.push(dst.len());

        if is_eol {
            break;
        }
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
        let mut record = Record::default();

        let mut src = &b"sq0\t0\t1\n"[..];
        read_record(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001");
        assert_eq!(record.0.bounds, Bounds::default());

        let mut src = &b"sq0\t0\t1\r\n"[..];
        read_record(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001");
        assert_eq!(record.0.bounds, Bounds::default());

        let mut src = &b"sq0\t0\t1\t.\n"[..];
        read_record(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001.");
        let mut bounds = Bounds::default();
        bounds.other_fields_ends.push(6);
        assert_eq!(record.0.bounds, bounds);

        Ok(())
    }
}
