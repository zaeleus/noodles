//! CSI index reader.

pub(crate) mod header;
pub mod reference_sequences;

use std::{
    error, fmt,
    io::{self, Read},
    num,
};

use byteorder::{LittleEndian, ReadBytesExt};

pub use self::header::read_header;
use self::{header::read_aux, reference_sequences::read_reference_sequences};
use crate::Index;

/// An error returned when a coordinate-sorted index fails to be read.
#[derive(Debug)]
pub enum ReadError {
    /// I/O error.
    Io(io::Error),
    /// The magic number is invalid.
    InvalidMagicNumber([u8; 4]),
    /// The min shift is invalid.
    InvalidMinShift(num::TryFromIntError),
    /// The depth is invalid.
    InvalidDepth(num::TryFromIntError),
    /// The header is invalid.
    InvalidHeader(header::ReadError),
    /// A reference sequence is invalid.
    InvalidReferenceSequences(reference_sequences::ReadError),
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidMagicNumber(_) => None,
            Self::InvalidMinShift(e) => Some(e),
            Self::InvalidDepth(e) => Some(e),
            Self::InvalidHeader(e) => Some(e),
            Self::InvalidReferenceSequences(e) => Some(e),
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidMagicNumber(_) => write!(f, "invalid magic number"),
            Self::InvalidMinShift(_) => write!(f, "invalid min shift"),
            Self::InvalidDepth(_) => write!(f, "invalid depth"),
            Self::InvalidHeader(_) => write!(f, "invalid header"),
            Self::InvalidReferenceSequences(_) => write!(f, "invalid reference sequences"),
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

pub(super) fn read_index<R>(reader: &mut R) -> Result<Index, ReadError>
where
    R: Read,
{
    read_magic(reader)?;

    let min_shift = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| u8::try_from(n).map_err(ReadError::InvalidMinShift))?;

    let depth = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| u8::try_from(n).map_err(ReadError::InvalidDepth))?;

    let header = read_aux(reader).map_err(ReadError::InvalidHeader)?;

    let reference_sequences =
        read_reference_sequences(reader, depth).map_err(ReadError::InvalidReferenceSequences)?;

    let n_no_coor = read_unplaced_unmapped_record_count(reader)?;

    let mut builder = Index::builder()
        .set_min_shift(min_shift)
        .set_depth(depth)
        .set_reference_sequences(reference_sequences);

    if let Some(hdr) = header {
        builder = builder.set_header(hdr);
    }

    if let Some(n_no_coor) = n_no_coor {
        builder = builder.set_unplaced_unmapped_record_count(n_no_coor);
    }

    Ok(builder.build())
}

fn read_magic<R>(reader: &mut R) -> Result<(), ReadError>
where
    R: Read,
{
    use crate::io::MAGIC_NUMBER;

    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(ReadError::InvalidMagicNumber(magic))
    }
}

fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> Result<Option<u64>, ReadError>
where
    R: Read,
{
    match reader.read_u64::<LittleEndian>() {
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(ReadError::Io(e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic() {
        let data = b"CSI\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());
    }

    #[test]
    fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ReadError::Io(e)) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"CSI";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ReadError::Io(e)) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ReadError::InvalidMagicNumber([b'M', b'T', b'h', b'd']))
        ));
    }

    #[test]
    fn test_read_unplaced_unmapped_record_count() -> Result<(), ReadError> {
        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, None);

        let data = [0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, Some(8));

        Ok(())
    }
}
