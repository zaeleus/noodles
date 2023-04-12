mod tag;
mod ty;
mod value;

pub(crate) use self::value::parse_value;

use std::io;

use self::{tag::parse_tag, ty::parse_type};
use crate::record::data::field::{Tag, Value};

pub(super) fn parse_field(src: &mut &[u8]) -> io::Result<(Tag, Value)> {
    use crate::reader::record::next_field;

    let mut buf = next_field(src);

    let tag = parse_tag(&mut buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    consume_delimiter(&mut buf)?;
    let ty = parse_type(&mut buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    consume_delimiter(&mut buf)?;
    let value =
        parse_value(&mut buf, ty).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((tag, value))
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
        let expected = (Tag::AlignmentHitCount, Value::from(1));
        assert_eq!(actual, expected);

        let actual = parse_field(&mut src)?;
        let expected = (Tag::Comment, Value::String(String::from("ndls")));
        assert_eq!(actual, expected);

        assert!(src.is_empty());

        Ok(())
    }
}
