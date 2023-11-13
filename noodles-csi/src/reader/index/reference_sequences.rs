mod bins;
mod metadata;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use self::{bins::read_bins, metadata::read_metadata};
use crate::index::ReferenceSequence;

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
    depth: u8,
) -> io::Result<Vec<ReferenceSequence>>
where
    R: Read,
{
    let n_ref = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    (0..n_ref)
        .map(|_| read_reference_sequence(reader, depth))
        .collect()
}

fn read_reference_sequence<R>(reader: &mut R, depth: u8) -> io::Result<ReferenceSequence>
where
    R: Read,
{
    let (bins, metadata) = read_bins(reader, depth)?;
    Ok(ReferenceSequence::new(bins, Vec::new(), metadata))
}
