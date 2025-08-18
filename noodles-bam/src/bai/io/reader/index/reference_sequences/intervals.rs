use std::io::{self, Read};

use noodles_bgzf as bgzf;

use crate::io::reader::num::{read_u32_le, read_u64_le};

pub(super) fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: Read,
{
    let n_intv = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioffset = read_u64_le(reader).map(bgzf::VirtualPosition::from)?;
        intervals.push(ioffset);
    }

    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_intervals() -> io::Result<()> {
        let src = [
            0x00, 0x00, 0x00, 0x00, // n_intv = 0
        ];
        assert!(read_intervals(&mut &src[..])?.is_empty());

        let src = [
            0x01, 0x00, 0x00, 0x00, // n_intv = 1
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[0] = 8
        ];
        assert_eq!(
            read_intervals(&mut &src[..])?,
            vec![bgzf::VirtualPosition::from(8)]
        );

        Ok(())
    }
}
