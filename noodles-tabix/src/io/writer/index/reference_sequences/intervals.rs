use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;

pub(super) fn write_intervals<W>(
    writer: &mut W,
    intervals: &[bgzf::VirtualPosition],
) -> io::Result<()>
where
    W: Write,
{
    let n_intv = i32::try_from(intervals.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_intv)?;

    for interval in intervals {
        let ioff = u64::from(*interval);
        writer.write_u64::<LittleEndian>(ioff)?;
    }

    Ok(())
}
