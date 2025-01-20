use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

pub(super) fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: Write,
{
    let n_chunk =
        u32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_chunk)?;

    for chunk in chunks {
        write_chunk(writer, chunk)?;
    }

    Ok(())
}

fn write_chunk<W>(writer: &mut W, chunk: &Chunk) -> io::Result<()>
where
    W: Write,
{
    let chunk_beg = u64::from(chunk.start());
    writer.write_u64::<LittleEndian>(chunk_beg)?;

    let chunk_end = u64::from(chunk.end());
    writer.write_u64::<LittleEndian>(chunk_end)?;

    Ok(())
}
