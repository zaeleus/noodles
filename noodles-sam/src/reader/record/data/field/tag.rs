use std::io;

use crate::record::data::field::Tag;

pub(super) fn parse_tag(src: &mut &[u8]) -> io::Result<Tag> {
    if src.len() < 2 {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let (buf, rest) = src.split_at(2);
    *src = rest;

    Tag::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() -> io::Result<()> {
        let mut src = &b"NH"[..];
        assert_eq!(parse_tag(&mut src)?, Tag::AlignmentHitCount);

        let mut src = &b""[..];
        assert!(matches!(
            parse_tag(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        let mut src = &b"X!"[..];
        assert!(matches!(
            parse_tag(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
