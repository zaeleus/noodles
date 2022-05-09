use std::io;

use bytes::Buf;
use noodles_sam::record::data::field::Tag;

pub fn get_tag<B>(src: &mut B) -> io::Result<Tag>
where
    B: Buf,
{
    let mut buf = [0; 2];

    if src.remaining() < buf.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    src.copy_to_slice(&mut buf);

    Tag::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_tag() -> io::Result<()> {
        let data = [b'N', b'H'];

        let mut reader = &data[..];
        let actual = get_tag(&mut reader)?;

        let expected = Tag::AlignmentHitCount;

        assert_eq!(actual, expected);

        Ok(())
    }
}
