mod chunks;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;

use self::chunks::write_chunks;
use super::write_metadata;
use crate::binning_index::index::reference_sequence::{index::BinnedIndex, Bin, Metadata};

pub(super) fn write_bins<W>(
    writer: &mut W,
    depth: u8,
    bins: &IndexMap<usize, Bin>,
    index: &BinnedIndex,
    metadata: Option<&Metadata>,
) -> io::Result<()>
where
    W: Write,
{
    let mut n_bin =
        i32::try_from(bins.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if metadata.is_some() {
        n_bin += 1;
    }

    writer.write_i32::<LittleEndian>(n_bin)?;

    for (id, bin) in bins {
        let first_record_start_position = index.get(id).copied().unwrap_or_default();
        write_bin(writer, *id, first_record_start_position, bin)?;
    }

    if let Some(m) = metadata {
        write_metadata(writer, depth, m)?;
    }

    Ok(())
}

fn write_bin<W>(
    writer: &mut W,
    id: usize,
    first_record_start_position: bgzf::VirtualPosition,
    bin: &Bin,
) -> io::Result<()>
where
    W: Write,
{
    let id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(id)?;

    let loffset = u64::from(first_record_start_position);
    writer.write_u64::<LittleEndian>(loffset)?;

    write_chunks(writer, bin.chunks())?;

    Ok(())
}
