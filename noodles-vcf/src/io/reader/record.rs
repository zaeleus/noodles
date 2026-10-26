use std::{
    io::{self, BufRead},
    mem,
};

use memchr::memchr2;

use crate::{Record, record::fields::Bounds};

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

const SEPARATOR: u8 = b'\t';

pub(crate) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    if is_eof(reader)? {
        return Ok(0);
    }

    let fields = record.fields_mut();

    let mut buf = mem::take(&mut fields.buf).into_bytes();
    buf.clear();

    let len = read_record_inner(reader, &mut buf, &mut fields.bounds)?;

    fields.buf =
        String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(len)
}

fn read_record_inner<R>(reader: &mut R, dst: &mut Vec<u8>, bounds: &mut Bounds) -> io::Result<usize>
where
    R: BufRead,
{
    let ends = [
        &mut bounds.reference_sequence_name_end,
        &mut bounds.variant_start_end,
        &mut bounds.ids_end,
        &mut bounds.reference_bases_end,
        &mut bounds.alternate_bases_end,
        &mut bounds.quality_score_end,
        &mut bounds.filters_end,
    ];

    let mut len = 0;

    for end in ends {
        len += read_required_field(reader, dst)?;
        *end = dst.len();
    }

    let (n, terminator) = read_last_required_field(reader, dst)?;
    len += n;
    bounds.info_end = dst.len();

    bounds.samples_keys_end = if matches!(terminator, Terminator::Separator) {
        len += read_required_field(reader, dst)?;

        let end = dst.len();

        len += match read_line(reader, dst)? {
            0 => return Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
            n => n,
        };

        if dst[end..].is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "missing fields"));
        }

        end
    } else {
        dst.len()
    };

    Ok(len)
}

fn is_eof<R>(reader: &mut R) -> io::Result<bool>
where
    R: BufRead,
{
    reader.fill_buf().map(|src| src.is_empty())
}

fn read_required_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    let prev_end = dst.len();
    let (len, terminator) = read_field(reader, dst)?;

    match terminator {
        Terminator::Separator => {
            if dst[prev_end..].is_empty() {
                Err(io::Error::new(io::ErrorKind::InvalidData, "missing field"))
            } else {
                Ok(len)
            }
        }
        Terminator::Eol => Err(io::Error::new(io::ErrorKind::InvalidData, "unexpected EOL")),
        Terminator::None => Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
    }
}

fn read_last_required_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<(usize, Terminator)>
where
    R: BufRead,
{
    let prev_end = dst.len();
    let (len, terminator) = read_field(reader, dst)?;

    if len == 0 {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    } else if dst[prev_end..].is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "missing field"));
    }

    Ok((len, terminator))
}

#[derive(Clone, Copy)]
enum Terminator {
    Separator,
    Eol,
    None,
}

impl Terminator {
    fn is_eol(self) -> bool {
        matches!(self, Terminator::Eol)
    }
}

fn read_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<(usize, Terminator)>
where
    R: BufRead,
{
    let mut terminator = Terminator::None;
    let mut len = 0;

    while matches!(terminator, Terminator::None) {
        let src = reader.fill_buf()?;

        if src.is_empty() {
            break;
        }

        let (buf, n) = match memchr2(SEPARATOR, LINE_FEED, src) {
            Some(i) => {
                terminator = match src[i] {
                    SEPARATOR => Terminator::Separator,
                    LINE_FEED => Terminator::Eol,
                    _ => unreachable!(),
                };

                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        dst.extend(buf);

        len += n;

        reader.consume(n);
    }

    if terminator.is_eol() && dst.ends_with(&[CARRIAGE_RETURN]) {
        dst.pop();
    }

    Ok((len, terminator))
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_until(LINE_FEED, buf)? {
        0 => Ok(0),
        n => {
            if buf.ends_with(&[LINE_FEED]) {
                buf.pop();

                if buf.ends_with(&[CARRIAGE_RETURN]) {
                    buf.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::fields::Bounds;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut record = Record::default();

        assert_eq!(read_record(&mut io::empty(), &mut record)?, 0);

        read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\n"[..], &mut record)?;
        assert_eq!(record, Record::default());

        read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\r\n"[..], &mut record)?;
        assert_eq!(record, Record::default());

        read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t."[..], &mut record)?;
        assert_eq!(record, Record::default());

        read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\tGT\t.\n"[..], &mut record)?;
        assert_eq!(record.fields().buf, ".1.N....GT.");
        assert_eq!(
            record.fields().bounds,
            Bounds {
                samples_keys_end: 10,
                ..Default::default()
            }
        );

        assert!(matches!(
            read_record(&mut &b"\n"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            read_record(&mut &b"\t"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            read_record(&mut &b".\t"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        assert!(matches!(
            read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        assert!(matches!(
            read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t\n"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\t"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\tGT"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\tGT\t"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_record(&mut &b".\t1\t.\tN\t.\t.\t.\t.\tGT\t\n"[..], &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_last_required_field() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        read_last_required_field(&mut &b"."[..], &mut buf)?;
        assert_eq!(buf, b".");

        buf.clear();
        assert!(matches!(
            read_last_required_field(&mut io::empty(), &mut buf),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        buf.clear();
        assert!(matches!(
            read_last_required_field(&mut &b"\t"[..], &mut buf),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        buf.clear();
        assert!(matches!(
            read_last_required_field(&mut &b"\n"[..], &mut buf),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_line() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, mut src: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_line(&mut src, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"noodles\n", b"noodles")?;
        t(&mut buf, b"noodles\r\n", b"noodles")?;
        t(&mut buf, b"noodles", b"noodles")?;

        Ok(())
    }
}
