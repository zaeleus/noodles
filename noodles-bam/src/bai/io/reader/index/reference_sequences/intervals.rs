use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

pub(super) fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: Read,
{
    let n_intv = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioffset = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioffset);
    }

    Ok(intervals)
}
