mod header;
mod magic_number;
mod reference_sequences;

use noodles_csi::BinningIndex;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::{
    header::write_header, magic_number::write_magic_number,
    reference_sequences::write_reference_sequences,
};
use crate::Index;

pub(super) async fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_magic_number(writer).await?;

    let reference_sequence_count = i32::try_from(index.reference_sequences().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(reference_sequence_count).await?;

    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing tabix header"))?;
    write_header(writer, header).await?;

    write_reference_sequences(writer, index.reference_sequences()).await?;

    if let Some(n) = index.unplaced_unmapped_record_count() {
        writer.write_u64_le(n).await?;
    }

    Ok(())
}
