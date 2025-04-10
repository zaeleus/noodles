mod chunks;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use indexmap::IndexMap;
use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};

use self::chunks::write_chunks;
use super::write_metadata;

pub(super) fn write_bins<W>(
    writer: &mut W,
    bins: &IndexMap<usize, Bin>,
    metadata: Option<&Metadata>,
) -> io::Result<()>
where
    W: Write,
{
    let mut n_bin =
        i32::try_from(bins.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if metadata.is_some() {
        n_bin = n_bin
            .checked_add(1)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "n_bin overflow"))?;
    }

    writer.write_i32::<LittleEndian>(n_bin)?;

    for (&id, bin) in bins {
        write_bin(writer, id, bin)?;
    }

    if let Some(metadata) = metadata {
        write_metadata(writer, metadata)?;
    }

    Ok(())
}

pub fn write_bin<W>(writer: &mut W, id: usize, bin: &Bin) -> io::Result<()>
where
    W: Write,
{
    let id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(id)?;
    write_chunks(writer, bin.chunks())?;
    Ok(())
}
