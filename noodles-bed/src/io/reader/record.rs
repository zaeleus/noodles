use std::io::{self, BufRead};

use memchr::{memchr, memchr2};

use crate::{Record, record::fields::Bounds};

pub(super) fn read_record_3<R>(reader: &mut R, record: &mut Record<3>) -> io::Result<usize>
where
    R: BufRead,
{
    let fields = &mut record.0;

    let dst = &mut fields.buf;
    dst.clear();

    let bounds = &mut fields.bounds;
    bounds.other_fields_ends.clear();

    skip_comment_lines(reader)?;

    let mut len = 0;

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[0] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[1] = dst.len();

    let (n, is_eol) = read_field(reader, dst)?;
    len += n;
    bounds.standard_fields_ends[2] = dst.len();

    if !is_eol {
        len += read_other_fields(reader, dst, bounds)?;
    }

    Ok(len)
}

pub(super) fn read_record_4<R>(reader: &mut R, record: &mut Record<4>) -> io::Result<usize>
where
    R: BufRead,
{
    let fields = &mut record.0;

    let dst = &mut fields.buf;
    dst.clear();

    let bounds = &mut fields.bounds;
    bounds.other_fields_ends.clear();

    skip_comment_lines(reader)?;

    let mut len = 0;

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[0] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[1] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[2] = dst.len();

    let (n, is_eol) = read_field(reader, dst)?;
    len += n;
    bounds.standard_fields_ends[3] = dst.len();

    if !is_eol {
        len += read_other_fields(reader, dst, bounds)?;
    }

    Ok(len)
}

pub(super) fn read_record_5<R>(reader: &mut R, record: &mut Record<5>) -> io::Result<usize>
where
    R: BufRead,
{
    let fields = &mut record.0;

    let dst = &mut fields.buf;
    dst.clear();

    let bounds = &mut fields.bounds;
    bounds.other_fields_ends.clear();

    skip_comment_lines(reader)?;

    let mut len = 0;

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[0] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[1] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[2] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[3] = dst.len();

    let (n, is_eol) = read_field(reader, dst)?;
    len += n;
    bounds.standard_fields_ends[4] = dst.len();

    if !is_eol {
        len += read_other_fields(reader, dst, bounds)?;
    }

    Ok(len)
}

pub(super) fn read_record_6<R>(reader: &mut R, record: &mut Record<6>) -> io::Result<usize>
where
    R: BufRead,
{
    let fields = &mut record.0;

    let dst = &mut fields.buf;
    dst.clear();

    let bounds = &mut fields.bounds;
    bounds.other_fields_ends.clear();

    skip_comment_lines(reader)?;

    let mut len = 0;

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[0] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[1] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[2] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[3] = dst.len();

    len += read_required_field(reader, dst)?;
    bounds.standard_fields_ends[4] = dst.len();

    let (n, is_eol) = read_field(reader, dst)?;
    len += n;
    bounds.standard_fields_ends[5] = dst.len();

    if !is_eol {
        len += read_other_fields(reader, dst, bounds)?;
    }

    Ok(len)
}

// § 1.3 "Terminology and concepts" (2022-01-05): "**comment line**: A **line** that starts with
// **#** with no horizontal whitespace beforehand."
const COMMENT_PREFIX: u8 = b'#';

fn skip_comment_lines<R>(reader: &mut R) -> io::Result<()>
where
    R: BufRead,
{
    loop {
        let src = reader.fill_buf()?;

        if src.starts_with(&[COMMENT_PREFIX]) {
            discard_line(reader)?;
        } else {
            break;
        }
    }

    Ok(())
}

fn discard_line<R>(reader: &mut R) -> io::Result<()>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';

    let mut is_eol = false;

    while !is_eol {
        let src = reader.fill_buf()?;

        if src.is_empty() {
            break;
        }

        let n = match memchr(LINE_FEED, src) {
            Some(i) => {
                is_eol = true;
                i + 1
            }
            None => src.len(),
        };

        reader.consume(n);
    }

    Ok(())
}

fn read_other_fields<R, const N: usize>(
    reader: &mut R,
    dst: &mut Vec<u8>,
    bounds: &mut Bounds<N>,
) -> io::Result<usize>
where
    R: BufRead,
{
    let mut len = 0;

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

        let (buf, n) = match memchr2(DELIMITER, LINE_FEED, src) {
            Some(i) => {
                r#match = Some(src[i]);
                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        dst.extend(buf);
        len += n;
        reader.consume(n);
    }

    let is_eol = matches!(r#match, Some(LINE_FEED));

    if is_eol && dst.ends_with(&[CARRIAGE_RETURN]) {
        dst.pop();
    }

    Ok((len, is_eol))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record_3() -> io::Result<()> {
        let mut record = Record::default();

        let mut src = &b"sq0\t0\t1\n"[..];
        read_record_3(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001");
        assert_eq!(record.0.bounds, Bounds::default());

        let mut src = &b"sq0\t0\t1\r\n"[..];
        read_record_3(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001");
        assert_eq!(record.0.bounds, Bounds::default());

        let mut src = &b"# noodles\nsq0\t0\t1\n"[..];
        read_record_3(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001");
        assert_eq!(record.0.bounds, Bounds::default());

        let mut src = &b"sq0\t0\t1\t.\n"[..];
        read_record_3(&mut src, &mut record)?;
        assert_eq!(record.0.buf, b"sq001.");
        let mut bounds = Bounds::default();
        bounds.other_fields_ends.push(6);
        assert_eq!(record.0.bounds, bounds);

        Ok(())
    }

    #[test]
    fn test_skip_comment_lines() -> io::Result<()> {
        let mut src = &b"sq0\t0\t1\n"[..];
        skip_comment_lines(&mut src)?;
        assert_eq!(src, &b"sq0\t0\t1\n"[..]);

        let mut src = &b"# noodles\nsq0\t0\t1\n"[..];
        skip_comment_lines(&mut src)?;
        assert_eq!(src, &b"sq0\t0\t1\n"[..]);

        let mut src = &b"# noodles\n# bed\nsq0\t0\t1\n"[..];
        skip_comment_lines(&mut src)?;
        assert_eq!(src, &b"sq0\t0\t1\n"[..]);

        Ok(())
    }

    #[test]
    fn test_discard_line() -> io::Result<()> {
        let mut src = &b"noo\ndles"[..];

        discard_line(&mut src)?;
        assert_eq!(src, &b"dles"[..]);

        discard_line(&mut src)?;
        assert!(src.is_empty());

        Ok(())
    }
}
