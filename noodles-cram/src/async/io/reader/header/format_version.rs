use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::file_definition::Version;

pub(super) async fn read_format_version<R>(reader: &mut R) -> io::Result<Version>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf).await?;
    Ok(Version::new(buf[0], buf[1]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_format_version() -> io::Result<()> {
        let src = [0x03, 0x00]; // (3, 0)
        let mut reader = &src[..];
        assert_eq!(read_format_version(&mut reader).await?, Version::new(3, 0));

        let src = [];
        let mut reader = &src[..];
        assert!(matches!(
            read_format_version(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
