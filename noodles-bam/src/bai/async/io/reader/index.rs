mod magic_number;
mod reference_sequences;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::{magic_number::read_magic_number, reference_sequences::read_reference_sequences};
use crate::bai::Index;

pub(super) async fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: AsyncRead + Unpin,
{
    read_magic_number(reader).await?;

    let reference_sequences = read_reference_sequences(reader).await?;
    let unplaced_unmapped_record_count = read_unplaced_unmapped_record_count(reader).await?;

    let mut builder = Index::builder().set_reference_sequences(reference_sequences);

    if let Some(n) = unplaced_unmapped_record_count {
        builder = builder.set_unplaced_unmapped_record_count(n);
    }

    Ok(builder.build())
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
    async fn test_read_unplaced_unmapped_record_count() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        assert_eq!(
            read_unplaced_unmapped_record_count(&mut reader).await?,
            None
        );

        let data = [0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(
            read_unplaced_unmapped_record_count(&mut reader).await?,
            Some(8)
        );

        Ok(())
    }
}
