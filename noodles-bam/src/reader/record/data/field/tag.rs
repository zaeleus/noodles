use std::{
    convert::TryFrom,
    io::{self, Read},
};

use noodles_sam::record::data::field::Tag;

pub fn read_tag<R>(reader: &mut R) -> io::Result<Tag>
where
    R: Read,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Tag::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_tag() -> io::Result<()> {
        let data = [b'N', b'H'];

        let mut reader = &data[..];
        let actual = read_tag(&mut reader)?;

        let expected = Tag::AlignmentHitCount;

        assert_eq!(actual, expected);

        Ok(())
    }
}
