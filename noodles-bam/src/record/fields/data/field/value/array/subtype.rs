use std::io;

use noodles_sam::alignment::record::data::field::value::array::Subtype;

pub(crate) fn decode_subtype(src: &mut &[u8]) -> io::Result<Subtype> {
    let Some((n, rest)) = src.split_first() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    match *n {
        b'c' => Ok(Subtype::Int8),
        b'C' => Ok(Subtype::UInt8),
        b's' => Ok(Subtype::Int16),
        b'S' => Ok(Subtype::UInt16),
        b'i' => Ok(Subtype::Int32),
        b'I' => Ok(Subtype::UInt32),
        b'f' => Ok(Subtype::Float),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid subtype",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_subtype() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Subtype) -> io::Result<()> {
            assert_eq!(decode_subtype(&mut src)?, expected);
            Ok(())
        }

        t(b"c", Subtype::Int8)?;
        t(b"C", Subtype::UInt8)?;
        t(b"s", Subtype::Int16)?;
        t(b"S", Subtype::UInt16)?;
        t(b"i", Subtype::Int32)?;
        t(b"I", Subtype::UInt32)?;
        t(b"f", Subtype::Float)?;

        let mut src = &[][..];
        assert!(matches!(
            decode_subtype(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"n"[..];
        assert!(matches!(
            decode_subtype(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
