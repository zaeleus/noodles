mod header;
mod magic_number;
mod reference_sequences;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::{
    header::read_header, magic_number::read_magic_number,
    reference_sequences::read_reference_sequences,
};
use crate::Index;

pub(super) async fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: AsyncRead + Unpin,
{
    read_magic_number(reader).await?;

    let reference_sequence_count = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let header = read_header(reader).await?;
    let reference_sequences = read_reference_sequences(reader, reference_sequence_count).await?;
    let unplaced_unmapped_record_count = read_unplaced_unmapped_record_count(reader).await?;

    let mut builder = Index::builder()
        .set_header(header)
        .set_reference_sequences(reference_sequences);

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
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}
