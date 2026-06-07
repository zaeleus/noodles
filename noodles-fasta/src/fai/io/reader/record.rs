use std::{
    io::{self, BufRead},
    num::NonZero,
};

use super::read_line;
use crate::fai::Record;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 5;

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
    if s.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty input"));
    }

    let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

    let name = parse_string(&mut fields)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing name"))?;

    let base_count = parse_u64(&mut fields)
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing base count"))?;

    let position = parse_u64(&mut fields)
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing position"))?;

    let line_base_count = parse_nonzero_u64(&mut fields)
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing line base count"))?;

    let line_width = parse_nonzero_u64(&mut fields)
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing line width"))?;

    Ok(Record::new(
        name,
        base_count,
        position,
        line_base_count,
        line_width,
    ))
}

fn parse_string<'a, I>(fields: &mut I) -> Option<&'a str>
where
    I: Iterator<Item = &'a str>,
{
    fields.next()
}

fn parse_u64<'a, I>(fields: &mut I) -> Option<io::Result<u64>>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().map(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn parse_nonzero_u64<'a, I>(fields: &mut I) -> Option<io::Result<NonZero<u64>>>
where
    I: Iterator<Item = &'a str>,
{
    fields.next().map(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_record() -> io::Result<()> {
        let line_base_count = const { NonZero::new(80).unwrap() };
        let line_width = const { NonZero::new(81).unwrap() };
        assert_eq!(
            parse_record("sq0\t10946\t4\t80\t81")?,
            Record::new("sq0", 10946, 4, line_base_count, line_width)
        );

        assert!(matches!(
            parse_record(""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_record("sq0"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_record("sq0\tndls"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
