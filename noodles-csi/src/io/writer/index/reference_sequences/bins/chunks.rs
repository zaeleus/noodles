use std::io::{self, Write};

use crate::{
    binning_index::index::reference_sequence::bin::Chunk,
    io::writer::num::{write_i32_le, write_u64_le},
};

pub(super) fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: Write,
{
    let n_chunk =
        i32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, n_chunk)?;

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

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;

    use super::*;

    #[test]
    fn test_write_chunks() -> io::Result<()> {
        let mut buf = Vec::new();
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
        );
        write_chunks(&mut buf, &[chunk])?;

        let expected = [
            0x01, 0x00, 0x00, 0x00, // n_chunk = 1
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_beg[0] = 8
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_end[0] = 13
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
