//! CSI index reader.

pub(crate) mod header;
mod magic_number;
pub mod reference_sequences;

use std::{
    error, fmt,
    io::{self, Read},
    num,
};

pub use self::header::read_header;
use self::{
    header::read_aux, magic_number::read_magic_number,
    reference_sequences::read_reference_sequences,
};
use super::{
    Index,
    num::{read_i32_le, read_u64_le},
};

/// An error returned when a coordinate-sorted index fails to be read.
#[derive(Debug)]
pub enum ReadError {
    /// I/O error.
    Io(io::Error),
    /// The magic number is invalid.
    InvalidMagicNumber(magic_number::ReadError),
    /// The min shift is invalid.
    InvalidMinShift(num::TryFromIntError),
    /// The depth is invalid.
    InvalidDepth,
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
            Self::InvalidDepth => None,
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
            Self::InvalidDepth => write!(f, "invalid depth"),
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
    read_magic_number(reader).map_err(ReadError::InvalidMagicNumber)?;

    let min_shift = read_min_shift(reader)?;
    let depth = read_depth(reader)?;

    let header = read_aux(reader).map_err(ReadError::InvalidHeader)?;

    let reference_sequences =
        read_reference_sequences(reader, depth).map_err(ReadError::InvalidReferenceSequences)?;

    let unplaced_unmapped_record_count = read_unplaced_unmapped_record_count(reader)?;

    let mut builder = Index::builder()
        .set_min_shift(min_shift)
        .set_depth(depth)
        .set_reference_sequences(reference_sequences);

    if let Some(hdr) = header {
        builder = builder.set_header(hdr);
    }

    if let Some(n) = unplaced_unmapped_record_count {
        builder = builder.set_unplaced_unmapped_record_count(n);
    }

    Ok(builder.build())
}

fn read_min_shift<R>(reader: &mut R) -> Result<u8, ReadError>
where
    R: Read,
{
    let n = read_i32_le(reader)?;
    u8::try_from(n).map_err(ReadError::InvalidMinShift)
}

fn read_depth<R>(reader: &mut R) -> Result<u8, ReadError>
where
    R: Read,
{
    const MAX_DEPTH: u8 = 9;

    let n = read_i32_le(reader)?;
    let depth = u8::try_from(n).map_err(|_| ReadError::InvalidDepth)?;

    if depth <= MAX_DEPTH {
        Ok(depth)
    } else {
        Err(ReadError::InvalidDepth)
    }
}

fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> Result<Option<u64>, ReadError>
where
    R: Read,
{
    match read_u64_le(reader) {
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(ReadError::Io(e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_min_shift() -> Result<(), ReadError> {
        let src = [0x0e, 0x00, 0x00, 0x00]; // min shift = 14
        assert_eq!(read_min_shift(&mut &src[..])?, 14);

        let src = [0xff, 0xff, 0xff, 0xff]; // min shift = -1
        assert!(matches!(
            read_min_shift(&mut &src[..]),
            Err(ReadError::InvalidMinShift(_))
        ));

        let src = [0x00, 0x01, 0x00, 0x00]; // min shift = 256
        assert!(matches!(
            read_min_shift(&mut &src[..]),
            Err(ReadError::InvalidMinShift(_))
        ));

        Ok(())
    }

    #[test]
    fn test_read_depth() -> Result<(), ReadError> {
        let src = [0x00, 0x00, 0x00, 0x00]; // depth = 0
        assert_eq!(read_depth(&mut &src[..])?, 0);

        let src = [0x05, 0x00, 0x00, 0x00]; // depth = 5
        assert_eq!(read_depth(&mut &src[..])?, 5);

        let src = [0x09, 0x00, 0x00, 0x00]; // depth = 9
        assert_eq!(read_depth(&mut &src[..])?, 9);

        let src = [0xff, 0xff, 0xff, 0xff]; // depth = -1
        assert!(matches!(
            read_depth(&mut &src[..]),
            Err(ReadError::InvalidDepth)
        ));

        let src = [0x0a, 0x00, 0x00, 0x00]; // depth = 10
        assert!(matches!(
            read_depth(&mut &src[..]),
            Err(ReadError::InvalidDepth)
        ));

        Ok(())
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
