mod subtype;
mod values;

use std::{
    io::{self, BufRead},
    mem,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::record::data::field::value::array::Subtype;

use self::subtype::decode_subtype;
pub use self::values::Values;

/// A raw BAM record data field array value.
#[derive(Debug, PartialEq)]
pub enum Array<'a> {
    /// An 8-bit integer array (`B:c`).
    Int8(Values<'a, i8>),
    /// An 8-bit unsigned integer array (`B:C`).
    UInt8(Values<'a, u8>),
    /// A 16-bit integer array (`B:s`).
    Int16(Values<'a, i16>),
    /// A 16-bit unsigned integer array (`B:S`).
    UInt16(Values<'a, u16>),
    /// A 32-bit integer array (`B:i`).
    Int32(Values<'a, i32>),
    /// A 32-bit unsigned integer array (`B:I`).
    UInt32(Values<'a, u32>),
    /// A single-precision floating-point array (`B:f`).
    Float(Values<'a, f32>),
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
            Ok(Array::Int8(Values::new(buf)))
        }
        Subtype::UInt8 => {
            let buf = &src[..n];
            src.consume(n);
            Ok(Array::UInt8(Values::new(buf)))
        }
        Subtype::Int16 => {
            let len = n * mem::size_of::<i16>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Int16(Values::new(buf)))
        }
        Subtype::UInt16 => {
            let len = n * mem::size_of::<u16>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::UInt16(Values::new(buf)))
        }
        Subtype::Int32 => {
            let len = n * mem::size_of::<i32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Int32(Values::new(buf)))
        }
        Subtype::UInt32 => {
            let len = n * mem::size_of::<u32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::UInt32(Values::new(buf)))
        }
        Subtype::Float => {
            let len = n * mem::size_of::<f32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Float(Values::new(buf)))
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

        t(
            &[b'c', 0x01, 0x00, 0x00, 0x00, 0x00],
            Array::Int8(Values::new(&[0x00])),
        )?;
        t(
            &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00],
            Array::UInt8(Values::new(&[0x00])),
        )?;
        t(
            &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::Int16(Values::new(&[0x00, 0x00])),
        )?;
        t(
            &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::UInt16(Values::new(&[0x00, 0x00])),
        )?;
        t(
            &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::Int32(Values::new(&[0x00, 0x00, 0x00, 0x00])),
        )?;
        t(
            &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::UInt32(Values::new(&[0x00, 0x00, 0x00, 0x00])),
        )?;
        t(
            &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Array::Float(Values::new(&[0x00, 0x00, 0x00, 0x00])),
        )?;

        Ok(())
    }
}
