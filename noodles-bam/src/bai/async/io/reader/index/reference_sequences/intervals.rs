use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: AsyncRead + Unpin,
{
    let n_intv = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioffset = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioffset);
    }

    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_intervals() -> io::Result<()> {
        let data = [
            0x03, 0x00, 0x00, 0x00, // n_intv = 3
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[0] = 8
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[1] = 13
            0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[2] = 21
        ];

        let mut reader = &data[..];
        let actual = read_intervals(&mut reader).await?;

        let expected = vec![
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
            bgzf::VirtualPosition::from(21),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
