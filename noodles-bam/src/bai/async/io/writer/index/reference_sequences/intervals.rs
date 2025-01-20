use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub(super) async fn write_intervals<W>(
    writer: &mut W,
    intervals: &[bgzf::VirtualPosition],
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_intv = u32::try_from(intervals.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_intv).await?;

    for &interval in intervals {
        let ioffset = u64::from(interval);
        writer.write_u64_le(ioffset).await?;
    }

    Ok(())
}
