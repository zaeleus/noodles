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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_chunks() -> io::Result<()> {
        let src = [
            0x01, 0x00, 0x00, 0x00, // n_chunk = 1
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_beg = 8
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_beg = 13
        ];
        let mut reader = &src[..];

        let actual = read_chunks(&mut reader)?;
        let expected = [Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
        )];

        assert_eq!(actual, expected);

        Ok(())
    }
}
