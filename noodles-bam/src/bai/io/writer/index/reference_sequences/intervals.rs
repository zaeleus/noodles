use std::io::{self, Write};

use noodles_bgzf as bgzf;

use crate::io::writer::num::{write_u32_le, write_u64_le};

pub(super) fn write_intervals<W>(
    writer: &mut W,
    intervals: &[bgzf::VirtualPosition],
) -> io::Result<()>
where
    W: Write,
{
    let n_intv = u32::try_from(intervals.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, n_intv)?;

    for interval in intervals {
        let ioffset = u64::from(*interval);
        write_u64_le(writer, ioffset)?;
    }

    Ok(())
}
