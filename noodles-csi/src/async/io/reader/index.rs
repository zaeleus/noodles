mod header;
mod magic_number;
mod reference_sequences;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::{
    header::read_aux, magic_number::read_magic_number,
    reference_sequences::read_reference_sequences,
};
use crate::{Index, io::MAX_DEPTH};

pub(super) async fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: AsyncRead + Unpin,
{
    read_magic_number(reader).await?;

    let min_shift = read_min_shift(reader).await?;

    let depth = reader
        .read_i32_le()
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    if depth > MAX_DEPTH {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid depth"));
    }

    let header = read_aux(reader).await?;

    let reference_sequences = read_reference_sequences(reader, depth).await?;
    let unplaced_unmapped_record_count = read_unplaced_unmapped_record_count(reader).await?;

    let mut builder = Index::builder()
        .set_min_shift(min_shift)
        .set_depth(depth)
        .set_reference_sequences(reference_sequences);

    if let Some(header) = header {
        builder = builder.set_header(header);
    }

    if let Some(n) = unplaced_unmapped_record_count {
        builder = builder.set_unplaced_unmapped_record_count(n);
    }

    Ok(builder.build())
}

async fn read_min_shift<R>(reader: &mut R) -> io::Result<u8>
where
    R: AsyncRead + Unpin,
{
    reader
        .read_i32_le()
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

async fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u64_le().await {
        Ok(n_no_coor) => Ok(Some(n_no_coor)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_min_shift() -> io::Result<()> {
        let src = [0x0e, 0x00, 0x00, 0x00]; // min shift = 14
        assert_eq!(read_min_shift(&mut &src[..]).await?, 14);

        let src = [0xff, 0xff, 0xff, 0xff]; // min shift = -1
        assert!(matches!(
            read_min_shift(&mut &src[..]).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let src = [0x00, 0x01, 0x00, 0x00]; // min shift = 256
        assert!(matches!(
            read_min_shift(&mut &src[..]).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
