mod subtype;

use std::io;

use self::subtype::parse_subtype;
use crate::record::data::field::value::array::Subtype;

#[derive(Debug, PartialEq)]
pub enum Array<'a> {
    Int8(&'a [u8]),
    UInt8(&'a [u8]),
    Int16(&'a [u8]),
    UInt16(&'a [u8]),
    Int32(&'a [u8]),
    UInt32(&'a [u8]),
    Float(&'a [u8]),
}

pub(super) fn parse_array<'a>(src: &mut &'a [u8]) -> io::Result<Array<'a>> {
    use super::parse_string;

    let subtype = parse_subtype(src)?;
    maybe_consume_delimiter(src)?;

    let buf = parse_string(src);

    match subtype {
        Subtype::Int8 => Ok(Array::Int8(buf)),
        Subtype::UInt8 => Ok(Array::UInt8(buf)),
        Subtype::Int16 => Ok(Array::Int16(buf)),
        Subtype::UInt16 => Ok(Array::UInt16(buf)),
        Subtype::Int32 => Ok(Array::Int32(buf)),
        Subtype::UInt32 => Ok(Array::UInt32(buf)),
        Subtype::Float => Ok(Array::Float(buf)),
    }
}

fn maybe_consume_delimiter(src: &mut &[u8]) -> io::Result<()> {
    const DELIMITER: u8 = b',';

    if let Some((b, rest)) = src.split_first() {
        if *b == DELIMITER {
            *src = rest;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid delimiter",
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_array() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Array<'_>) -> io::Result<()> {
            assert_eq!(parse_array(&mut src)?, expected);
            Ok(())
        }

        t(b"c,0", Array::Int8(b"0"))?;
        t(b"C,0", Array::UInt8(b"0"))?;
        t(b"s,0", Array::Int16(b"0"))?;
        t(b"S,0", Array::UInt16(b"0"))?;
        t(b"i,0", Array::Int32(b"0"))?;
        t(b"I,0", Array::UInt32(b"0"))?;
        t(b"f,0.0", Array::Float(b"0.0"))?;

        Ok(())
    }
}
