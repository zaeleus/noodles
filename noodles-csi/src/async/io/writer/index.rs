mod header;
mod magic_number;
mod reference_sequences;

use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::{
    header::write_header, magic_number::write_magic_number,
    reference_sequences::write_reference_sequences,
};
use crate::{BinningIndex, Index};

pub(super) async fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_magic_number(writer).await?;

    write_header(writer, index.min_shift(), index.depth(), index.header()).await?;
    write_reference_sequences(writer, index.depth(), index.reference_sequences()).await?;

    if let Some(n) = index.unplaced_unmapped_record_count() {
        writer.write_u64_le(n).await?;
    }

    Ok(())
}
