use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_file_id<R>(reader: &mut R) -> io::Result<[u8; 20]>
where
    R: AsyncRead + Unpin,
{
    let mut file_id = [0; 20];
    reader.read_exact(&mut file_id).await?;
    Ok(file_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_file_id() -> io::Result<()> {
        let src = [
            0x00, 0xac, 0x24, 0xf8, 0xc4, 0x2d, 0xc2, 0xa5, 0x56, 0xa0, 0x85, 0x1c, 0xa5, 0xef,
            0xf0, 0xfc, 0x6d, 0x40, 0x33, 0x4d,
        ];
        let mut reader = &src[..];
        assert_eq!(read_file_id(&mut reader).await?, src);

        let src = [];
        let mut reader = &src[..];
        assert!(matches!(
            read_file_id(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}