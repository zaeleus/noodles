#![doc(hidden)]

pub mod bins;
mod metadata;

use std::{
    error, fmt,
    io::{self, Read},
    num,
};

use byteorder::{LittleEndian, ReadBytesExt};

use self::bins::read_bins;
pub use self::metadata::read_metadata;
use crate::index::{reference_sequence::index::BinnedIndex, ReferenceSequence};

/// An error returned when CSI reference sequences fail to be read.
#[derive(Debug)]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The reference sequence count is invalid.
    InvalidReferenceSequenceCount(num::TryFromIntError),
    /// The bins are invalid.
    InvalidBins(bins::ReadError),
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidReferenceSequenceCount(e) => Some(e),
            Self::InvalidBins(e) => Some(e),
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidReferenceSequenceCount(_) => write!(f, "invalid reference sequence count"),
            Self::InvalidBins(_) => write!(f, "invalid bins"),
        }
    }
}

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
    depth: u8,
) -> Result<Vec<ReferenceSequence<BinnedIndex>>, ReadError>
where
    R: Read,
{
    let n_ref = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| usize::try_from(n).map_err(ReadError::InvalidReferenceSequenceCount))?;

    (0..n_ref)
        .map(|_| read_reference_sequence(reader, depth))
        .collect()
}

fn read_reference_sequence<R>(
    reader: &mut R,
    depth: u8,
) -> Result<ReferenceSequence<BinnedIndex>, ReadError>
where
    R: Read,
{
    let (bins, index, metadata) = read_bins(reader, depth).map_err(ReadError::InvalidBins)?;
    Ok(ReferenceSequence::new(bins, index, metadata))
}
