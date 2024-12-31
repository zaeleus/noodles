use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_format_version<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf).await?;
    Ok((buf[0], buf[1]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_format_version() -> io::Result<()> {
        let data = [0x02, 0x01];
        let mut reader = &data[..];
        assert_eq!(read_format_version(&mut reader).await?, (2, 1));
        Ok(())
    }
}
