use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::bai::MAGIC_NUMBER;

pub(super) async fn read_magic_number<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let mut magic = [0; MAGIC_NUMBER.len()];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAI header",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic_number() {
        let data = b"BAI\x01";
        let mut reader = &data[..];
        assert!(read_magic_number(&mut reader).await.is_ok());
    }

    #[tokio::test]
    async fn test_read_magic_number_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAI";
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
