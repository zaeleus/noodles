use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: AsyncRead + Unpin,
{
    // n_intv
    let interval_count = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(interval_count);

    for _ in 0..interval_count {
        // ioffset
        let offset = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(offset);
    }

    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_intervals() -> io::Result<()> {
        let src = [
            0x00, 0x00, 0x00, 0x00, // n_intv = 0
        ];
        assert!(read_intervals(&mut &src[..]).await?.is_empty());

        let src = [
            0x01, 0x00, 0x00, 0x00, // n_intv = 1
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[0] = 8
        ];
        assert_eq!(
            read_intervals(&mut &src[..]).await?,
            vec![bgzf::VirtualPosition::from(8)]
        );

        Ok(())
    }
}
