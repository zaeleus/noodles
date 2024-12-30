use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(crate) async fn read_format_version<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: AsyncRead + Unpin,
{
    let major_version = reader.read_u8().await?;
    let minor_version = reader.read_u8().await?;

    Ok((major_version, minor_version))
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
