mod subtype;
mod values;

use std::io;

use self::subtype::parse_subtype;
pub use self::values::Values;
use crate::alignment::record::data::field::value::{array::Subtype, Array};

pub(super) fn parse_array<'a>(src: &mut &'a [u8]) -> io::Result<Array<'a>> {
    use super::parse_string;

    let subtype = parse_subtype(src)?;
    maybe_consume_delimiter(src)?;

    let buf = parse_string(src);

    match subtype {
        Subtype::Int8 => Ok(Array::Int8(Box::new(Values::new(buf)))),
        Subtype::UInt8 => Ok(Array::UInt8(Box::new(Values::new(buf)))),
        Subtype::Int16 => Ok(Array::Int16(Box::new(Values::new(buf)))),
        Subtype::UInt16 => Ok(Array::UInt16(Box::new(Values::new(buf)))),
        Subtype::Int32 => Ok(Array::Int32(Box::new(Values::new(buf)))),
        Subtype::UInt32 => Ok(Array::UInt32(Box::new(Values::new(buf)))),
        Subtype::Float => Ok(Array::Float(Box::new(Values::new(buf)))),
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
        let mut src = &b"c,0"[..];
        if let Array::Int8(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &b"C,0"[..];
        if let Array::UInt8(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &b"s,0"[..];
        if let Array::Int16(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &b"S,0"[..];
        if let Array::UInt16(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &b"i,0"[..];
        if let Array::Int32(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &b"I,0"[..];
        if let Array::UInt32(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &b"f,0.0"[..];
        if let Array::Float(values) = parse_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0.0]);
        } else {
            panic!();
        }

        Ok(())
    }
}
