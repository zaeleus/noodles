use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: AsyncRead + Unpin,
{
    let n_intv = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioff = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioff);
    }

    Ok(intervals)
}
