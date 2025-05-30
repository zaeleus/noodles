use std::{
    error, fmt,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use crate::binning_index::index::reference_sequence::{Metadata, bin::METADATA_CHUNK_COUNT};

/// An error returned when CSI reference sequence metadata fail to be read.
#[derive(Debug)]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The chunk count is invalid.
    InvalidChunkCount(u32),
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ReadError::Io(e) => Some(e),
            ReadError::InvalidChunkCount(_) => None,
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadError::Io(_) => write!(f, "I/O error"),
            ReadError::InvalidChunkCount(actual) => write!(
                f,
                "invalid chunk count: expected {METADATA_CHUNK_COUNT}, got {actual}"
            ),
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

pub fn read_metadata<R>(reader: &mut R) -> Result<Metadata, ReadError>
where
    R: Read,
{
    let n_chunk = reader.read_u32::<LittleEndian>()?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(ReadError::InvalidChunkCount(n_chunk));
    }

    let ref_beg = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let ref_end = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let n_mapped = reader.read_u64::<LittleEndian>()?;
    let n_unmapped = reader.read_u64::<LittleEndian>()?;

    Ok(Metadata::new(ref_beg, ref_end, n_mapped, n_unmapped))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_metadata() -> Result<(), ReadError> {
        let data = [
            0x02, 0x00, 0x00, 0x00, // n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        let mut reader = &data[..];
        let actual = read_metadata(&mut reader)?;

        let expected = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_metadata_with_invalid_chunk_count() {
        let data = [
            0x01, 0x00, 0x00, 0x00, // n_chunk = 1
        ];
        let mut reader = &data[..];

        assert!(matches!(
            read_metadata(&mut reader),
            Err(ReadError::InvalidChunkCount(1))
        ));
    }
}
