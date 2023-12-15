mod subtype;
mod values;

use std::{
    io::{self, BufRead},
    mem,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::{
    alignment::record::data::field::value::Array, record::data::field::value::array::Subtype,
};

use self::subtype::decode_subtype;
pub use self::values::Values;

pub(super) fn decode_array<'a>(src: &mut &'a [u8]) -> io::Result<Array<'a>> {
    let subtype = decode_subtype(src)?;

    let n = src.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    match subtype {
        Subtype::Int8 => {
            let buf = &src[..n];
            src.consume(n);
            Ok(Array::Int8(Box::new(Values::new(buf))))
        }
        Subtype::UInt8 => {
            let buf = &src[..n];
            src.consume(n);
            Ok(Array::UInt8(Box::new(Values::new(buf))))
        }
        Subtype::Int16 => {
            let len = n * mem::size_of::<i16>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Int16(Box::new(Values::new(buf))))
        }
        Subtype::UInt16 => {
            let len = n * mem::size_of::<u16>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::UInt16(Box::new(Values::new(buf))))
        }
        Subtype::Int32 => {
            let len = n * mem::size_of::<i32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Int32(Box::new(Values::new(buf))))
        }
        Subtype::UInt32 => {
            let len = n * mem::size_of::<u32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::UInt32(Box::new(Values::new(buf))))
        }
        Subtype::Float => {
            let len = n * mem::size_of::<f32>();
            let buf = &src[..len];
            src.consume(len);
            Ok(Array::Float(Box::new(Values::new(buf))))
        }
    }
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
