use std::io;

const LENGTH: usize = 2;

/// A raw BAM record data field tag.
pub type Tag = [u8; LENGTH];

pub(super) fn parse_tag(src: &mut &[u8]) -> io::Result<Tag> {
    if src.len() < LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let (buf, rest) = src.split_at(LENGTH);

    // SAFETY: `buf` is 2 bytes.
    let tag = buf.try_into().unwrap();

    *src = rest;

    Ok(tag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() -> io::Result<()> {
        let mut src = &b"NH"[..];
        assert_eq!(parse_tag(&mut src)?, [b'N', b'H']);
        Ok(())
    }
}
