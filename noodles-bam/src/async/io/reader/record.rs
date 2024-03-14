use std::mem;

use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_record<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let block_size = match read_block_size(reader).await? {
        0 => return Ok(0),
        n => n,
    };

    buf.resize(block_size, 0);
    reader.read_exact(buf).await?;

    Ok(block_size)
}

async fn read_block_size<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; mem::size_of::<u32>()];
    read_exact_or_eof(reader, &mut buf).await?;
    let n = u32::from_le_bytes(buf);
    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_exact_or_eof<R>(reader: &mut R, mut buf: &mut [u8]) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let mut bytes_read = 0;

    while !buf.is_empty() {
        match reader.read(buf).await {
            Ok(0) => break,
            Ok(n) => {
                buf = &mut buf[n..];
                bytes_read += n;
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    if bytes_read > 0 && !buf.is_empty() {
        Err(io::Error::new(io::ErrorKind::UnexpectedEof, "early eof"))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_block_size() -> io::Result<()> {
        let data = [0x08, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader).await?, 8);

        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader).await?, 0);

        let data = [0x08];
        let mut reader = &data[..];
        assert!(matches!(
            read_block_size(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
