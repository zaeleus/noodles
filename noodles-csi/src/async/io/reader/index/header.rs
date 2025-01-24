use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::binning_index::index::Header;

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<(u8, u8, Option<Header>)>
where
    R: AsyncRead + Unpin,
{
    let min_shift = reader
        .read_i32_le()
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let depth = reader
        .read_i32_le()
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let header = read_aux(reader).await?;

    Ok((min_shift, depth, header))
}

async fn read_aux<R>(reader: &mut R) -> io::Result<Option<Header>>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::index::header::read_header as read_tabix_header;

    let l_aux = reader.read_i32_le().await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if l_aux > 0 {
        let mut aux = vec![0; l_aux];
        reader.read_exact(&mut aux).await?;

        let mut rdr = &aux[..];
        read_tabix_header(&mut rdr)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Some)
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_header() -> io::Result<()> {
        let data = [
            0x0e, 0x00, 0x00, 0x00, // min_shift = 14
            0x05, 0x00, 0x00, 0x00, // depth = 5
            0x00, 0x00, 0x00, 0x00, // l_aux = 0
        ];

        let mut reader = &data[..];
        let (min_shift, depth, header) = read_header(&mut reader).await?;

        assert_eq!(min_shift, 14);
        assert_eq!(depth, 5);
        assert!(header.is_none());

        Ok(())
    }
}
