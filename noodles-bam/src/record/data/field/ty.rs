use std::io;

use noodles_sam::alignment::record::data::field::Type;

pub(crate) fn decode_type(src: &mut &[u8]) -> io::Result<Type> {
    let Some((n, rest)) = src.split_first() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    match n {
        b'A' => Ok(Type::Character),
        b'c' => Ok(Type::Int8),
        b'C' => Ok(Type::UInt8),
        b's' => Ok(Type::Int16),
        b'S' => Ok(Type::UInt16),
        b'i' => Ok(Type::Int32),
        b'I' => Ok(Type::UInt32),
        b'f' => Ok(Type::Float),
        b'Z' => Ok(Type::String),
        b'H' => Ok(Type::Hex),
        b'B' => Ok(Type::Array),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_type() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Type) -> io::Result<()> {
            assert_eq!(decode_type(&mut src)?, expected);
            Ok(())
        }

        t(b"A", Type::Character)?;
        t(b"c", Type::Int8)?;
        t(b"C", Type::UInt8)?;
        t(b"s", Type::Int16)?;
        t(b"S", Type::UInt16)?;
        t(b"i", Type::Int32)?;
        t(b"I", Type::UInt32)?;
        t(b"f", Type::Float)?;
        t(b"Z", Type::String)?;
        t(b"H", Type::Hex)?;
        t(b"B", Type::Array)?;

        let mut src = &[][..];
        assert!(matches!(
            decode_type(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"n"[..];
        assert!(matches!(
            decode_type(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
