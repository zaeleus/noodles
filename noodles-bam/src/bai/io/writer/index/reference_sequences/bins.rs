mod chunks;

use std::io::{self, Write};

use indexmap::IndexMap;
use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};

use self::chunks::write_chunks;
use super::write_metadata;
use crate::io::writer::num::write_u32_le;

pub(super) fn write_bins<W>(
    writer: &mut W,
    bins: &IndexMap<usize, Bin>,
    metadata: Option<&Metadata>,
) -> io::Result<()>
where
    W: Write,
{
    let n_bin = u32::try_from(bins.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .and_then(|n| {
            if metadata.is_some() {
                n.checked_add(1)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "n_bin overflow"))
            } else {
                Ok(n)
            }
        })?;

    write_u32_le(writer, n_bin)?;

    for (&id, bin) in bins {
        write_bin(writer, id, bin)?;
    }

    if let Some(metadata) = metadata {
        write_metadata(writer, metadata)?;
    }

    Ok(())
}

fn write_bin<W>(writer: &mut W, id: usize, bin: &Bin) -> io::Result<()>
where
    W: Write,
{
    let id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, id)?;
    write_chunks(writer, bin.chunks())?;
    Ok(())
}
