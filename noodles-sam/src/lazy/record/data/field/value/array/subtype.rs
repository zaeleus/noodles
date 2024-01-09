use std::io;

use crate::alignment::record::data::field::value::array::Subtype;

pub(super) fn parse_subtype(src: &mut &[u8]) -> io::Result<Subtype> {
    let (b, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    let subtype = match *b {
        b'c' => Subtype::Int8,
        b'C' => Subtype::UInt8,
        b's' => Subtype::Int16,
        b'S' => Subtype::UInt16,
        b'i' => Subtype::Int32,
        b'I' => Subtype::UInt32,
        b'f' => Subtype::Float,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid subtype",
            ));
        }
    };

    *src = rest;

    Ok(subtype)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_subtype() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Subtype) -> io::Result<()> {
            assert_eq!(parse_subtype(&mut src)?, expected);
            Ok(())
        }

        t(b"c", Subtype::Int8)?;
        t(b"C", Subtype::UInt8)?;
        t(b"s", Subtype::Int16)?;
        t(b"S", Subtype::UInt16)?;
        t(b"i", Subtype::Int32)?;
        t(b"I", Subtype::UInt32)?;
        t(b"f", Subtype::Float)?;

        let mut src = &b""[..];
        assert!(matches!(
            parse_subtype(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"n"[..];
        assert!(matches!(
            parse_subtype(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
