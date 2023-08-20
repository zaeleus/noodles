use std::io::{self, Read};

use crate::lazy::record::data::Tag;

pub(super) fn decode_tag(src: &mut &[u8]) -> io::Result<Tag> {
    let mut buf = [0; 2];
    src.read_exact(&mut buf)?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_tag() -> io::Result<()> {
        let mut src = &[b'N', b'H'][..];
        assert_eq!(decode_tag(&mut src)?, [b'N', b'H']);
        Ok(())
    }
}
