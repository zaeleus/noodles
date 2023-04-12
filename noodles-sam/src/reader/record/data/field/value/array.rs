use std::io;

use super::subtype::parse_subtype;
use crate::record::data::field::value::{array::Subtype, Array};

pub(super) fn parse_array(src: &mut &[u8]) -> io::Result<Array> {
    const DELIMITER: u8 = b',';

    fn consume_delimiter(src: &mut &[u8]) -> io::Result<()> {
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

    let subtype = parse_subtype(src)?;

    match subtype {
        Subtype::Int8 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::Int8(values))
        }
        Subtype::UInt8 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::UInt8(values))
        }
        Subtype::Int16 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::Int16(values))
        }
        Subtype::UInt16 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::UInt16(values))
        }
        Subtype::Int32 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::Int32(values))
        }
        Subtype::UInt32 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::UInt32(values))
        }
        Subtype::Float => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Array::Float(values))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_array() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Array) -> io::Result<()> {
            let actual = parse_array(&mut src)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(b"c", Array::Int8(vec![]))?;
        t(b"c,0", Array::Int8(vec![0]))?;
        t(b"c,0,0", Array::Int8(vec![0, 0]))?;

        t(b"C", Array::UInt8(vec![]))?;
        t(b"C,0", Array::UInt8(vec![0]))?;
        t(b"C,0,0", Array::UInt8(vec![0, 0]))?;

        t(b"s", Array::Int16(vec![]))?;
        t(b"s,0", Array::Int16(vec![0]))?;
        t(b"s,0,0", Array::Int16(vec![0, 0]))?;

        t(b"S", Array::UInt16(vec![]))?;
        t(b"S,0", Array::UInt16(vec![0]))?;
        t(b"S,0,0", Array::UInt16(vec![0, 0]))?;

        t(b"i", Array::Int32(vec![]))?;
        t(b"i,0", Array::Int32(vec![0]))?;
        t(b"i,0,0", Array::Int32(vec![0, 0]))?;

        t(b"I", Array::UInt32(vec![]))?;
        t(b"I,0", Array::UInt32(vec![0]))?;
        t(b"I,0,0", Array::UInt32(vec![0, 0]))?;

        t(b"f", Array::Float(vec![]))?;
        t(b"f,0", Array::Float(vec![0.0]))?;
        t(b"f,0,0", Array::Float(vec![0.0, 0.0]))?;

        Ok(())
    }
}
