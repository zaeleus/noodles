use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_magic_number<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    use crate::MAGIC_NUMBER;

    let mut magic = [0; 4];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid tabix header",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic_number() {
        let data = b"TBI\x01";
        let mut reader = &data[..];
        assert!(read_magic_number(&mut reader).await.is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
