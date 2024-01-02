use std::io::{self, Read};

/// A raw BAM record data field tag.
pub type Tag = [u8; 2];

pub(crate) fn decode_tag(src: &mut &[u8]) -> io::Result<Tag> {
    let mut buf = [0; 2];
    src.read_exact(&mut buf)?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_tag() -> io::Result<()> {
        let mut src = &[b'N', b'H'][..];
        assert_eq!(decode_tag(&mut src)?, [b'N', b'H']);
        Ok(())
    }
}
