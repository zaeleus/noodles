mod tag;
mod value;

pub(crate) use self::value::parse_value;

use std::io;

use self::tag::parse_tag;
use crate::record::data::Field;

pub(super) fn parse_field(src: &mut &[u8]) -> io::Result<Option<Field>> {
    use crate::reader::record::next_field;

    let mut buf = next_field(src);

    let tag = match parse_tag(&mut buf) {
        Ok(t) => t,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    };

    consume_delimiter(&mut buf)?;
    let ty = value::parse_type(&mut buf)?;

    consume_delimiter(&mut buf)?;
    let value = parse_value(&mut buf, ty)?;

    Ok(Some(Field::new(tag, value)))
}

fn consume_delimiter(src: &mut &[u8]) -> io::Result<()> {
    const DELIMITER: u8 = b':';

    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    if *n == DELIMITER {
        Ok(())
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        use crate::record::data::field::{Tag, Value};

        let mut src = &b"NH:i:1\tCO:Z:ndls"[..];

        let actual = parse_field(&mut src)?;
        let expected = Field::new(Tag::AlignmentHitCount, Value::from(1u8));
        assert_eq!(actual, Some(expected));

        let actual = parse_field(&mut src)?;
        let expected = Field::new(Tag::Comment, Value::String(String::from("ndls")));
        assert_eq!(actual, Some(expected));

        assert!(parse_field(&mut src)?.is_none());

        Ok(())
    }
}
