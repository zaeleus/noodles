mod subtype;

use std::{
    io::{self, BufRead},
    mem,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::record::data::field::value::array::Subtype;

use self::subtype::decode_subtype;

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

pub(super) fn decode_array<'a>(src: &mut &'a [u8]) -> io::Result<Array<'a>> {
    let subtype = decode_subtype(src)?;

    let n = src.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    match subtype {
        Subtype::Int8 => {
            let buf = &src[..n];
            src.consume(n);
            Ok(Array::Int8(buf))
        }
        Subtype::UInt8 => {
            let buf = &src[..n];
            src.consume(n);
            Ok(Array::UInt8(buf))
        }
        Subtype::Int16 => {
            let len = n * mem::size_of::<i16>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Int16(buf))
        }
        Subtype::UInt16 => {
            let len = n * mem::size_of::<u16>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::UInt16(buf))
        }
        Subtype::Int32 => {
            let len = n * mem::size_of::<i32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Int32(buf))
        }
        Subtype::UInt32 => {
            let len = n * mem::size_of::<u32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::UInt32(buf))
        }
        Subtype::Float => {
            let len = n * mem::size_of::<f32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Float(buf))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_array() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Array<'_>) -> io::Result<()> {
            assert_eq!(decode_array(&mut src)?, expected);
            Ok(())
        }

        t(&[b'c', 0x01, 0x00, 0x00, 0x00, 0x00], Array::Int8(&[0x00]))?;
        t(&[b'C', 0x01, 0x00, 0x00, 0x00, 0x00], Array::UInt8(&[0x00]))?;
        t(
            &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::Int16(&[0x00, 0x00]),
        )?;
        t(
            &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::UInt16(&[0x00, 0x00]),
        )?;
        t(
            &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::Int32(&[0x00, 0x00, 0x00, 0x00]),
        )?;
        t(
            &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::UInt32(&[0x00, 0x00, 0x00, 0x00]),
        )?;
        t(
            &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::Float(&[0x00, 0x00, 0x00, 0x00]),
        )?;

        Ok(())
    }
}
