use std::io;

use crate::alignment::record::data::field::Tag;

pub(super) fn parse_tag(src: &mut &[u8]) -> io::Result<Tag> {
    let ([a, b], rest) =
        split_first_chunk(src).ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;
    let tag = Tag::new(*a, *b);
    *src = rest;
    Ok(tag)
}

// TODO: Use `slice::split_first_chunk` when the MSRV is raised to or above Rust 1.77.0.
fn split_first_chunk<const N: usize>(src: &[u8]) -> Option<(&[u8; N], &[u8])> {
    if src.len() < N {
        None
    } else {
        // SAFETY: `src.len` >= `N`.
        let (head, tail) = src.split_at(N);
        <&[u8; N]>::try_from(head).ok().map(|chunk| (chunk, tail))
    }
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
