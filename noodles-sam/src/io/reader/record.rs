use std::io::{self, BufRead};

use super::read_line;
use crate::Record;

use memchr::memchr2;

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

    let buf = &mut fields.buf;
    buf.clear();

    let bounds = &mut fields.bounds;

    let ends = [
        &mut bounds.name_end,
        &mut bounds.flags_end,
        &mut bounds.reference_sequence_name_end,
        &mut bounds.alignment_start_end,
        &mut bounds.mapping_quality_end,
        &mut bounds.cigar_end,
        &mut bounds.mate_reference_sequence_name_end,
        &mut bounds.mate_alignment_start_end,
        &mut bounds.template_length_end,
        &mut bounds.sequence_end,
    ];

    let mut len = 0;

    for end in ends {
        len += read_required_field(reader, buf)?;
        *end = buf.len();
    }

    let (n, terminator) = read_last_required_field(reader, buf)?;
    len += n;
    bounds.quality_scores_end = buf.len();

    if matches!(terminator, Terminator::Separator) {
        len += read_line(reader, buf)?;
    }

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
    let (len, terminator) = read_field(reader, dst)?;

    match terminator {
        Terminator::Separator => Ok(len),
        Terminator::Eol => Err(io::Error::new(io::ErrorKind::InvalidData, "unexpected EOL")),
        Terminator::None => Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
    }
}

fn read_last_required_field<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<(usize, Terminator)>
where
    R: BufRead,
{
    read_field(reader, dst)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::fields::Bounds;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let mut record = Record::default();

        assert_eq!(read_record(&mut io::empty(), &mut record)?, 0);

        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"[..];
        read_record(&mut src, &mut record)?;
        assert_eq!(record.fields().buf, b"*4*0255**00**");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\r\n"[..];
        read_record(&mut src, &mut record)?;
        assert_eq!(record.fields().buf, b"*4*0255**00**");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tNH:i:0\n"[..];
        read_record(&mut src, &mut record)?;
        assert_eq!(record.fields().buf, b"*4*0255**00**NH:i:0");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"\n"[..];
        assert!(matches!(
            read_record(&mut src, &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
