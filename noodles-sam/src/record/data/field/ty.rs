use std::io;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Type {
    Character,
    Integer,
    Float,
    String,
    Hex,
    Array,
}

pub(super) fn parse_type(src: &mut &[u8]) -> io::Result<Type> {
    if let Some((b, rest)) = src.split_first() {
        let ty = match b {
            b'A' => Type::Character,
            b'i' => Type::Integer,
            b'f' => Type::Float,
            b'Z' => Type::String,
            b'H' => Type::Hex,
            b'B' => Type::Array,
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
        };

        *src = rest;

        Ok(ty)
    } else {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_type() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Type) -> io::Result<()> {
            assert_eq!(parse_type(&mut src)?, expected);
            Ok(())
        }

        t(b"A", Type::Character)?;
        t(b"i", Type::Integer)?;
        t(b"f", Type::Float)?;
        t(b"Z", Type::String)?;
        t(b"H", Type::Hex)?;
        t(b"B", Type::Array)?;

        let mut src = &b""[..];
        assert!(matches!(
            parse_type(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &b"n"[..];
        assert!(matches!(
            parse_type(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
