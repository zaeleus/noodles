mod magic_number;
mod reference_sequences;

use noodles_csi::BinningIndex;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::{magic_number::write_magic_number, reference_sequences::write_reference_sequences};
use crate::bai::Index;

pub(super) async fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_magic_number(writer).await?;
    write_reference_sequences(writer, index.reference_sequences()).await?;

    if let Some(unplaced_unmapped_record_count) = index.unplaced_unmapped_record_count() {
        write_unplaced_unmapped_record_count(writer, unplaced_unmapped_record_count).await?;
    }

    Ok(())
}

async fn write_unplaced_unmapped_record_count<W>(
    writer: &mut W,
    unplaced_unmapped_record_count: u64,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_u64_le(unplaced_unmapped_record_count).await
}
