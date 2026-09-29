use std::{
    error, fmt,
    io::{self, Read},
};

use noodles_bgzf as bgzf;

use crate::{
    binning_index::index::reference_sequence::{Metadata, bin::METADATA_CHUNK_COUNT},
    io::reader::num::{read_u32_le, read_u64_le},
};

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
    read_chunk_count(reader)?;

    let ref_beg = read_u64_le(reader).map(bgzf::VirtualPosition::from)?;
    let ref_end = read_u64_le(reader).map(bgzf::VirtualPosition::from)?;

    let n_mapped = read_u64_le(reader)?;
    let n_unmapped = read_u64_le(reader)?;

    Ok(Metadata::new(ref_beg, ref_end, n_mapped, n_unmapped))
}

fn read_chunk_count<R>(reader: &mut R) -> Result<(), ReadError>
where
    R: Read,
{
    let n = read_u32_le(reader)?;

    if n == METADATA_CHUNK_COUNT {
        Ok(())
    } else {
        Err(ReadError::InvalidChunkCount(n))
    }
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
    fn test_read_chunk_count() {
        assert!(read_chunk_count(&mut &[0x02, 0x00, 0x00, 0x00][..]).is_ok());

        assert!(matches!(
            read_chunk_count(&mut io::empty()),
            Err(ReadError::Io(e)) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
        assert!(matches!(
            read_chunk_count(&mut &[0x00, 0x00, 0x00, 0x00][..]),
            Err(ReadError::InvalidChunkCount(0))
        ));
        assert!(matches!(
            read_chunk_count(&mut &[0x08, 0x00, 0x00, 0x00][..]),
            Err(ReadError::InvalidChunkCount(8))
        ));
    }
}
