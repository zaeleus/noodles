use noodles_sam::record::data::field::Tag;

use tokio::io::{self, AsyncRead, AsyncReadExt};

pub async fn read_tag<R>(reader: &mut R) -> io::Result<Tag>
where
    R: AsyncRead + Unpin,
{
    use std::str;

    let mut buf = [0; 2];
    reader.read_exact(&mut buf).await?;

    str::from_utf8(&buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_tag() -> io::Result<()> {
        let data = [b'N', b'H'];

        let mut reader = &data[..];
        let actual = read_tag(&mut reader).await?;

        let expected = Tag::AlignmentHitCount;

        assert_eq!(actual, expected);

        Ok(())
    }
}
