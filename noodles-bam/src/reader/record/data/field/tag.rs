use std::io;

use bytes::Buf;
use noodles_sam::alignment::record::data::field::Tag;

pub fn get_tag<B>(src: &mut B) -> io::Result<Tag>
where
    B: Buf,
{
    if src.remaining() < 2 {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let buf = [src.get_u8(), src.get_u8()];
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
