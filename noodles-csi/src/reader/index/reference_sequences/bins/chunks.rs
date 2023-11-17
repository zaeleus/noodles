use std::{
    error, fmt,
    io::{self, Read},
    num,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use crate::index::reference_sequence::bin::Chunk;

/// An error returned when CSI reference sequence bin chunks fail to be read.
#[derive(Debug)]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The chunk count is invalid.
    InvalidChunkCount(num::TryFromIntError),
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidChunkCount(e) => Some(e),
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidChunkCount(_) => write!(f, "invalid chunk count"),
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

pub fn read_chunks<R>(reader: &mut R) -> Result<Vec<Chunk>, ReadError>
where
    R: Read,
{
    let n_chunk = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| usize::try_from(n).map_err(ReadError::InvalidChunkCount))?;

    (0..n_chunk).map(|_| read_chunk(reader)).collect()
}

fn read_chunk<R>(reader: &mut R) -> Result<Chunk, ReadError>
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
    fn test_read_chunks() -> Result<(), ReadError> {
        let src = [
            0x01, 0x00, 0x00, 0x00, // n_chunk = 1
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_beg[0] = 8
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_end[0] = 13
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
