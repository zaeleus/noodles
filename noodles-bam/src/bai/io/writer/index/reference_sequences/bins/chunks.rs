use std::io::{self, Write};

use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::io::writer::num::{write_u32_le, write_u64_le};

pub(super) fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: Write,
{
    let n_chunk =
        u32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, n_chunk)?;

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
    write_u64_le(writer, chunk_beg)?;

    let chunk_end = u64::from(chunk.end());
    write_u64_le(writer, chunk_end)?;

    Ok(())
}
