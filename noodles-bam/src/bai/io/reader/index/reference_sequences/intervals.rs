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
