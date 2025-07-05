use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::io::MAGIC_NUMBER;

pub(super) async fn read_magic_number<R>(reader: &mut R) -> io::Result<[u8; MAGIC_NUMBER.len()]>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; MAGIC_NUMBER.len()];
    reader.read_exact(&mut buf).await?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic_number() -> io::Result<()> {
        let src = b"BAM\x01";
        let mut reader = &src[..];
        assert_eq!(read_magic_number(&mut reader).await?, *b"BAM\x01");

        let mut src = &[][..];
        assert!(matches!(
            read_magic_number(&mut src).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
