mod tag;
mod ty;
pub mod value;

use std::io;

use self::ty::Type;
use self::{tag::parse_tag, ty::parse_type, value::parse_value};
use crate::alignment::record::data::field::{Tag, Value};

pub(super) fn parse_field<'a>(src: &mut &'a [u8]) -> io::Result<(Tag, Value<'a>)> {
    let tag = parse_tag(src)?;
    consume_delimiter(src)?;
    let ty = parse_type(src)?;
    consume_delimiter(src)?;
    let value = parse_value(src, ty)?;
    maybe_consume_terminator(src)?;

    Ok((tag, value))
}

fn consume_delimiter(src: &mut &[u8]) -> io::Result<()> {
    const DELIMITER: u8 = b':';

    if let Some((b, rest)) = src.split_first() {
        if *b == DELIMITER {
            *src = rest;
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid delimiter",
            ))
        }
    } else {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    }
}

fn maybe_consume_terminator(src: &mut &[u8]) -> io::Result<()> {
    const TERMINATOR: u8 = b'\t';

    if let Some((b, rest)) = src.split_first() {
        if *b == TERMINATOR {
            *src = rest;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid field terminator",
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        let mut src = &b"NH:i:1"[..];
        assert!(matches!(
            parse_field(&mut src)?,
            (Tag::ALIGNMENT_HIT_COUNT, Value::Int32(1))
        ));

        let mut src = &b"NH:i:1\t"[..];
        assert!(matches!(
            parse_field(&mut src)?,
            (Tag::ALIGNMENT_HIT_COUNT, Value::Int32(1))
        ));

        let mut src = &b""[..];
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"NH_i:1"[..];
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let mut src = &b"NH:i_1"[..];
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let mut src = &b"NH:i:1\n"[..];
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
