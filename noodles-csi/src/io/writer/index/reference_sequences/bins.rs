mod chunks;

use std::io::{self, Write};

use indexmap::IndexMap;
use noodles_bgzf as bgzf;

use self::chunks::write_chunks;
use super::write_metadata;
use crate::{
    binning_index::index::reference_sequence::{Bin, Metadata, index::BinnedIndex, parent_id},
    io::writer::num::{write_i32_le, write_u32_le, write_u64_le},
};

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

    write_i32_le(writer, n_bin)?;

    for (&id, bin) in bins {
        let first_record_start_position = first_record_start_position(index, id);
        write_bin(writer, id, first_record_start_position, bin)?;
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
    write_u32_le(writer, id)?;

    let loffset = u64::from(first_record_start_position);
    write_u64_le(writer, loffset)?;

    write_chunks(writer, bin.chunks())?;

    Ok(())
}

fn first_record_start_position(index: &BinnedIndex, mut id: usize) -> bgzf::VirtualPosition {
    let mut min_position = index.get(&id).copied().unwrap_or_default();

    while let Some(pid) = parent_id(id)
        && let Some(position) = index.get(&pid)
    {
        if *position < min_position {
            min_position = *position;
        }

        id = pid;
    }

    min_position
}
