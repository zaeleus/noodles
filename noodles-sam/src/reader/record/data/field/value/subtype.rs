use std::io;

use crate::record::data::field::value::Subtype;

pub(super) fn parse_subtype(src: &mut &[u8]) -> io::Result<Subtype> {
    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Subtype::try_from(*n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_subtype() -> io::Result<()> {
        let mut src = &b"i"[..];
        assert_eq!(parse_subtype(&mut src)?, Subtype::Int32);

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
