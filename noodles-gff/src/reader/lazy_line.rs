use std::{
    io::{self, BufRead},
    mem, str,
};

use super::read_line;
use crate::lazy;

pub(crate) fn read_lazy_line<R>(reader: &mut R, line: &mut lazy::Line) -> io::Result<usize>
where
    R: BufRead,
{
    const DEFAULT_LINE: lazy::Line = lazy::Line::Comment(String::new());
    const DIRECTIVE_PREFIX: &str = "##";

    let prev_line = mem::replace(line, DEFAULT_LINE);
    let mut buf: String = prev_line.into();

    match peek_line_type(reader)? {
        Some(LineType::Comment) => {
            buf.clear();

            let n = read_line(reader, &mut buf)?;

            *line = if buf.starts_with(DIRECTIVE_PREFIX) {
                lazy::Line::Directive(buf)
            } else {
                lazy::Line::Comment(buf)
            };

            Ok(n)
        }
        Some(LineType::Record) => {
            let (n, bounds) = read_lazy_record(reader, &mut buf)?;
            *line = lazy::Line::Record(lazy::Record { buf, bounds });
            Ok(n)
        }
        None => Ok(0),
    }
}

enum LineType {
    Comment,
    Record,
}

fn peek_line_type<R>(reader: &mut R) -> io::Result<Option<LineType>>
where
    R: BufRead,
{
    const COMMENT_PREFIX: u8 = b'#';

    let src = reader.fill_buf()?;

    Ok(src.first().map(|&b| match b {
        COMMENT_PREFIX => LineType::Comment,
        _ => LineType::Record,
    }))
}

fn read_lazy_record<R>(
    reader: &mut R,
    buf: &mut String,
) -> io::Result<(usize, lazy::record::Bounds)>
where
    R: BufRead,
{
    buf.clear();

    let mut len = 0;
    let mut bounds = lazy::record::Bounds::default();

    len += read_field(reader, buf)?;
    bounds.reference_sequence_name_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.source_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.type_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.start_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.end_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.score_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.strand_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.phase_end = buf.len();

    len += read_line(reader, buf)?;

    Ok((len, bounds))
}

fn read_field<R>(reader: &mut R, dst: &mut String) -> io::Result<usize>
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

        let (buf, n) = match src.iter().position(|&b| b == DELIMITER) {
            Some(i) => {
                is_delimiter = true;
                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        let s = str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        dst.push_str(s);

        len += n;

        reader.consume(n);
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lazy::record::Bounds;

    #[test]
    fn test_read_lazy_line() -> io::Result<()> {
        let mut line = lazy::Line::default();

        let mut src = &b"##gff-version 3"[..];
        read_lazy_line(&mut src, &mut line)?;
        assert_eq!(line, lazy::Line::Directive(String::from("##gff-version 3")));

        let mut src = &b"#noodles"[..];
        read_lazy_line(&mut src, &mut line)?;
        assert_eq!(line, lazy::Line::Comment(String::from("#noodles")));

        let mut src = &b".\t.\t.\t1\t1\t.\t.\t.\t."[..];
        read_lazy_line(&mut src, &mut line)?;
        assert_eq!(
            line,
            lazy::Line::Record(lazy::Record {
                buf: String::from("...11...."),
                bounds: Bounds::default()
            })
        );

        Ok(())
    }
}
