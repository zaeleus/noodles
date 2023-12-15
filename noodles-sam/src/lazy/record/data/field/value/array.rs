mod subtype;
mod values;

use std::io;

use self::subtype::parse_subtype;
pub use self::values::Values;
use crate::record::data::field::value::array::Subtype;

#[derive(Debug, PartialEq)]
pub enum Array<'a> {
    Int8(Values<'a, i8>),
    UInt8(Values<'a, u8>),
    Int16(Values<'a, i16>),
    UInt16(Values<'a, u16>),
    Int32(Values<'a, i32>),
    UInt32(Values<'a, u32>),
    Float(Values<'a, f32>),
}

pub(super) fn parse_array<'a>(src: &mut &'a [u8]) -> io::Result<Array<'a>> {
    use super::parse_string;

    let subtype = parse_subtype(src)?;
    maybe_consume_delimiter(src)?;

    let buf = parse_string(src);

    match subtype {
        Subtype::Int8 => Ok(Array::Int8(Values::new(buf))),
        Subtype::UInt8 => Ok(Array::UInt8(Values::new(buf))),
        Subtype::Int16 => Ok(Array::Int16(Values::new(buf))),
        Subtype::UInt16 => Ok(Array::UInt16(Values::new(buf))),
        Subtype::Int32 => Ok(Array::Int32(Values::new(buf))),
        Subtype::UInt32 => Ok(Array::UInt32(Values::new(buf))),
        Subtype::Float => Ok(Array::Float(Values::new(buf))),
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

        t(b"c,0", Array::Int8(Values::new(b"0")))?;
        t(b"C,0", Array::UInt8(Values::new(b"0")))?;
        t(b"s,0", Array::Int16(Values::new(b"0")))?;
        t(b"S,0", Array::UInt16(Values::new(b"0")))?;
        t(b"i,0", Array::Int32(Values::new(b"0")))?;
        t(b"I,0", Array::UInt32(Values::new(b"0")))?;
        t(b"f,0.0", Array::Float(Values::new(b"0.0")))?;

        Ok(())
    }
}
