use std::mem;

use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::Record;

pub(super) async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let block_size = u32::try_from(record.block_size())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(block_size).await?;

    writer.write_i32_le(record.ref_id).await?;
    writer.write_i32_le(record.pos).await?;

    let l_read_name = u8::try_from(record.read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u8(l_read_name).await?;

    let mapq = u8::from(record.mapping_quality());
    writer.write_u8(mapq).await?;

    writer.write_u16_le(record.bin()).await?;

    let n_cigar_op = u16::try_from(record.cigar().len() / mem::size_of::<u32>())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16_le(n_cigar_op).await?;

    let flag = u16::from(record.flags());
    writer.write_u16_le(flag).await?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_seq).await?;

    writer.write_i32_le(record.next_ref_id).await?;
    writer.write_i32_le(record.next_pos).await?;

    writer.write_i32_le(record.template_length()).await?;

    writer.write_all(&record.read_name).await?;

    for &raw_op in record.cigar().as_ref().iter() {
        writer.write_u32_le(raw_op).await?;
    }

    writer.write_all(record.sequence().as_ref()).await?;
    writer.write_all(record.quality_scores().as_ref()).await?;
    writer.write_all(record.data().as_ref()).await?;

    Ok(())
}
