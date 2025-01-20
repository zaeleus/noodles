use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::bai::MAGIC_NUMBER;

pub(super) async fn write_magic_number<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(&MAGIC_NUMBER).await
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_magic_number() -> io::Result<()> {
        let mut buf = Vec::new();
        write_magic_number(&mut buf).await?;
        assert_eq!(buf, b"BAI\x01");
        Ok(())
    }
}
