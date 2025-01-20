use std::io::{self, Write};

use noodles_bgzf as bgzf;

use byteorder::{LittleEndian, WriteBytesExt};

pub(super) fn write_intervals<W>(
    writer: &mut W,
    intervals: &[bgzf::VirtualPosition],
) -> io::Result<()>
where
    W: Write,
{
    let n_intv = u32::try_from(intervals.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_intv)?;

    for interval in intervals {
        let ioffset = u64::from(*interval);
        writer.write_u64::<LittleEndian>(ioffset)?;
    }

    Ok(())
}
