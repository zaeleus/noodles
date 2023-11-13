use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use crate::index::reference_sequence::bin::Chunk;

pub(super) fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    (0..n_chunk).map(|_| read_chunk(reader)).collect()
}

fn read_chunk<R>(reader: &mut R) -> io::Result<Chunk>
where
    R: Read,
{
    let chunk_beg = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let chunk_end = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    Ok(Chunk::new(chunk_beg, chunk_end))
}
