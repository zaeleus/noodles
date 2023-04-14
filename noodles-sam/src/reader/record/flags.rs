use std::io;

use crate::record::Flags;

pub(crate) fn parse_flags(src: &[u8]) -> io::Result<Flags> {
    lexical_core::parse::<u16>(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(Flags::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_flags() -> io::Result<()> {
        assert_eq!(parse_flags(b"0")?, Flags::empty());
        assert_eq!(parse_flags(b"4")?, Flags::UNMAPPED);

        assert!(matches!(
            parse_flags(b""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_flags(b"-4"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_flags(b"n"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
