use std::io;

use crate::alignment::record::data::field::Tag;

pub(super) fn parse_tag(src: &mut &[u8]) -> io::Result<Tag> {
    let ([b0, b1], rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(Tag::new(*b0, *b1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() -> io::Result<()> {
        let mut src = &b"NH"[..];
        assert_eq!(parse_tag(&mut src)?, Tag::ALIGNMENT_HIT_COUNT);
        Ok(())
    }
}
