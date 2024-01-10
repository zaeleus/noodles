use std::io;

use crate::alignment::record::data::field::Tag;

pub(super) fn parse_tag(src: &mut &[u8]) -> io::Result<Tag> {
    if src.len() < 2 {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let (buf, rest) = src.split_at(2);

    // SAFETY: `buf` is 2 bytes.
    let raw_tag: [u8; 2] = buf.try_into().unwrap();
    let tag = Tag::from(raw_tag);

    *src = rest;

    Ok(tag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() -> io::Result<()> {
        use crate::alignment::record::data::field::tag;

        let mut src = &b"NH"[..];
        assert_eq!(parse_tag(&mut src)?, tag::ALIGNMENT_HIT_COUNT);
        Ok(())
    }
}
