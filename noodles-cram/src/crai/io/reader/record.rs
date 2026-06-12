use std::io::{self, BufRead};

use noodles_core::Position;

use super::read_line;
use crate::crai::Record;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 6;

pub(super) fn read_record<R>(
    reader: &mut R,
    buf: &mut String,
    record: &mut Record,
) -> io::Result<usize>
where
    R: BufRead,
{
    match read_line(reader, buf)? {
        0 => Ok(0),
        n => {
            *record = parse_record(buf)?;
            Ok(n)
        }
    }
}

pub(crate) fn parse_record(s: &str) -> io::Result<Record> {
    const UNMAPPED: i32 = -1;

    let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

    let reference_sequence_id = parse_i32(&mut fields).and_then(|n| match n {
        UNMAPPED => Ok(None),
        _ => usize::try_from(n)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    })?;

    let alignment_start = parse_position(&mut fields)?;
    let alignment_span = parse_span(&mut fields)?;
    let offset = parse_u64(&mut fields)?;
    let landmark = parse_u64(&mut fields)?;
    let slice_length = parse_u64(&mut fields)?;

    Ok(Record::new(
        reference_sequence_id,
        alignment_start,
        alignment_span,
        offset,
        landmark,
        slice_length,
    ))
}

fn parse_i32<'a, I>(fields: &mut I) -> io::Result<i32>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

fn parse_u64<'a, I>(fields: &mut I) -> io::Result<u64>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

fn parse_position<'a, I>(fields: &mut I) -> io::Result<Option<Position>>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .map(Position::new)
}

fn parse_span<'a, I>(fields: &mut I) -> io::Result<usize>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_record() -> io::Result<()> {
        let actual = parse_record("0\t10946\t6765\t17711\t233\t317811")?;
        let expected = Record::new(Some(0), Position::new(10946), 6765, 17711, 233, 317811);
        assert_eq!(actual, expected);

        assert!(matches!(
            parse_record("0\t10946"),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            parse_record("0\t10946\tnoodles"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_record("-8\t10946\t6765\t17711\t233\t317811"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
