mod subtype;
mod values;

use std::{io, mem};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::{
    alignment::record::data::field::value::Array, record::data::field::value::array::Subtype,
};

use self::subtype::decode_subtype;
pub use self::values::Values;

pub(super) fn decode_array<'a>(src: &mut &'a [u8]) -> io::Result<Array<'a>> {
    let subtype = decode_subtype(src)?;
    let buf = decode_raw_array(src, subtype)?;

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

fn decode_length(src: &mut &[u8]) -> io::Result<usize> {
    src.read_u32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

fn decode_raw_array<'a>(src: &mut &'a [u8], subtype: Subtype) -> io::Result<&'a [u8]> {
    let n = decode_length(src)?;

    let len = match subtype {
        Subtype::Int8 => n,
        Subtype::UInt8 => n,
        Subtype::Int16 => n * mem::size_of::<i16>(),
        Subtype::UInt16 => n * mem::size_of::<u16>(),
        Subtype::Int32 => n * mem::size_of::<i32>(),
        Subtype::UInt32 => n * mem::size_of::<u32>(),
        Subtype::Float => n * mem::size_of::<f32>(),
    };

    let (buf, rest) = src.split_at(len);

    *src = rest;

    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_array() -> io::Result<()> {
        let mut src = &[b'c', 0x01, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::Int8(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::UInt8(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::Int16(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::UInt16(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::Int32(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::UInt32(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        let mut src = &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00][..];
        if let Array::Float(values) = decode_array(&mut src)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0.0]);
        } else {
            panic!();
        }

        Ok(())
    }
}
