use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::binning_index::index::Header;

pub(super) async fn read_aux<R>(reader: &mut R) -> io::Result<Option<Header>>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::index::header::read_header as read_tabix_header;

    let len = reader.read_i32_le().await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if len > 0 {
        let mut aux = vec![0; len];
        reader.read_exact(&mut aux).await?;

        let mut aux_reader = &aux[..];
        read_tabix_header(&mut aux_reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Some)
    } else {
        Ok(None)
    }
}
