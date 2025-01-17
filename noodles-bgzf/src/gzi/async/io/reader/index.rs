use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::gzi::Index;

pub(super) async fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: AsyncRead + Unpin,
{
    let len = reader.read_u64_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut offsets = Vec::with_capacity(len);

    for _ in 0..len {
        let compressed = reader.read_u64_le().await?;
        let uncompressed = reader.read_u64_le().await?;
        offsets.push((compressed, uncompressed));
    }

    match reader.read_u8().await {
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected trailing data",
        )),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(Index::from(offsets)),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_index() -> io::Result<()> {
        let src = [
            0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // len = 2
            0x3c, 0x12, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 4668
            0x2e, 0x53, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 21294
            0x02, 0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // compressed_offset = 23810
            0x01, 0x52, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // uncompressed_offset = 86529
        ];

        let mut reader = &src[..];

        assert_eq!(
            read_index(&mut reader).await?,
            Index::from(vec![(4668, 21294), (23810, 86529)])
        );

        Ok(())
    }

    #[tokio::test]
    async fn test_read_index_with_no_entries() -> io::Result<()> {
        let src = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // len = 0
        let mut reader = &src[..];
        assert_eq!(read_index(&mut reader).await?, Index::default());
        Ok(())
    }

    #[tokio::test]
    async fn test_read_index_with_fewer_than_len_entries() -> io::Result<()> {
        let src = [0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // len = 1
        let mut reader = &src[..];

        assert!(matches!(
            read_index(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_index_with_trailing_data() -> io::Result<()> {
        let src = [
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // len = 0
            0x00,
        ];

        let mut reader = &src[..];

        assert!(matches!(
            read_index(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
