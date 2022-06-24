use std::io;

use crate::record::data::field::value::Type;

pub(crate) fn parse_type(src: &mut &[u8]) -> io::Result<Type> {
    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Type::try_from(*n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_type() -> io::Result<()> {
        let mut src = &b"i"[..];
        assert_eq!(parse_type(&mut src)?, Type::Int32);

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
