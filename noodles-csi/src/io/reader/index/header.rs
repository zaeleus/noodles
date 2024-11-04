mod reference_sequence_names;

use std::{
    error, fmt,
    io::{self, Read},
    num,
};

use byteorder::{LittleEndian, ReadBytesExt};

use self::reference_sequence_names::read_reference_sequence_names;
use crate::binning_index::index::{header::format, Header};

/// An error returned when a CSI header fails to be read.
#[derive(Debug)]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The aux length is invalid.
    InvalidAuxLength(num::TryFromIntError),
    /// The header format is invalid.
    InvalidFormat(format::TryFromIntError),
    /// The header reference sequence index is invalid.
    InvalidReferenceSequenceNameIndex(num::TryFromIntError),
    /// The header reference sequence index value is invalid.
    InvalidReferenceSequenceNameIndexValue,
    /// The header start position index is invalid.
    InvalidStartPositionIndex(num::TryFromIntError),
    /// The header start position index value is invalid.
    InvalidStartPositionIndexValue,
    /// The header end position index is invalid.
    InvalidEndPositionIndex(num::TryFromIntError),
    /// The header line comment prefix is invalid.
    InvalidLineCommentPrefix(num::TryFromIntError),
    /// The header line skip count is invalid.
    InvalidLineSkipCount(num::TryFromIntError),
    /// The reference sequence names are invalid.
    InvalidReferenceSequenceNames(reference_sequence_names::ReadError),
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidAuxLength(e) => Some(e),
            Self::InvalidFormat(e) => Some(e),
            Self::InvalidReferenceSequenceNameIndex(e) => Some(e),
            Self::InvalidStartPositionIndex(e) => Some(e),
            Self::InvalidEndPositionIndex(e) => Some(e),
            Self::InvalidLineCommentPrefix(e) => Some(e),
            Self::InvalidLineSkipCount(e) => Some(e),
            Self::InvalidReferenceSequenceNames(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidAuxLength(_) => write!(f, "invalid aux length"),
            Self::InvalidFormat(_) => write!(f, "invalid format"),
            Self::InvalidReferenceSequenceNameIndex(_) => {
                write!(f, "invalid reference sequence name index")
            }
            Self::InvalidReferenceSequenceNameIndexValue => {
                write!(f, "invalid reference sequence name index value")
            }
            Self::InvalidStartPositionIndex(_) => write!(f, "invalid start position index"),
            Self::InvalidStartPositionIndexValue => write!(f, "invalid start position index value"),
            Self::InvalidEndPositionIndex(_) => write!(f, "invalid end position index"),
            Self::InvalidLineCommentPrefix(_) => write!(f, "invalid line comment prefix"),
            Self::InvalidLineSkipCount(_) => write!(f, "invalid line skip count"),
            Self::InvalidReferenceSequenceNames(_) => write!(f, "invalid reference sequence names"),
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

pub(super) fn read_aux<R>(reader: &mut R) -> Result<Option<Header>, ReadError>
where
    R: Read,
{
    let l_aux = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| u64::try_from(n).map_err(ReadError::InvalidAuxLength))?;

    if l_aux > 0 {
        let mut aux_reader = reader.take(l_aux);
        read_header(&mut aux_reader).map(Some)
    } else {
        Ok(None)
    }
}

#[doc(hidden)]
pub fn read_header<R>(reader: &mut R) -> Result<Header, ReadError>
where
    R: Read,
{
    use crate::binning_index::index::header::Format;

    let format =
        read_i32_le(reader).and_then(|n| Format::try_from(n).map_err(ReadError::InvalidFormat))?;

    let col_seq = read_reference_sequence_name_index(reader)?;
    let col_beg = read_start_position_index(reader)?;
    let col_end = read_end_position_index(reader)?;

    let meta = read_i32_le(reader)
        .and_then(|b| u8::try_from(b).map_err(ReadError::InvalidLineCommentPrefix))?;

    let skip = read_i32_le(reader)
        .and_then(|n| u32::try_from(n).map_err(ReadError::InvalidLineSkipCount))?;

    let names =
        read_reference_sequence_names(reader).map_err(ReadError::InvalidReferenceSequenceNames)?;

    Ok(Header::builder()
        .set_format(format)
        .set_reference_sequence_name_index(col_seq)
        .set_start_position_index(col_beg)
        .set_end_position_index(col_end)
        .set_line_comment_prefix(meta)
        .set_line_skip_count(skip)
        .set_reference_sequence_names(names)
        .build())
}

fn read_i32_le<R>(reader: &mut R) -> Result<i32, ReadError>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>().map_err(ReadError::Io)
}

fn read_reference_sequence_name_index<R>(reader: &mut R) -> Result<usize, ReadError>
where
    R: Read,
{
    read_i32_le(reader).and_then(|i| {
        usize::try_from(i)
            .map_err(ReadError::InvalidReferenceSequenceNameIndex)
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or(ReadError::InvalidReferenceSequenceNameIndexValue)
            })
    })
}

fn read_start_position_index<R>(reader: &mut R) -> Result<usize, ReadError>
where
    R: Read,
{
    read_i32_le(reader).and_then(|i| {
        usize::try_from(i)
            .map_err(ReadError::InvalidStartPositionIndex)
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or(ReadError::InvalidStartPositionIndexValue)
            })
    })
}

fn read_end_position_index<R>(reader: &mut R) -> Result<Option<usize>, ReadError>
where
    R: Read,
{
    read_i32_le(reader).and_then(|i| match i {
        0 => Ok(None),
        _ => usize::try_from(i)
            .map(|n| {
                // SAFETY: `n` is > 0.
                n - 1
            })
            .map(Some)
            .map_err(ReadError::InvalidEndPositionIndex),
    })
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

    #[test]
    fn test_read_aux() -> Result<(), ReadError> {
        let src = [0x00, 0x00, 0x00, 0x00];
        let mut reader = &src[..];
        assert!(read_aux(&mut reader)?.is_none());

        let src = [
            0x24, 0x00, 0x00, 0x00, // l_aux = 36
            0x02, 0x00, 0x00, 0x00, // format = 2 (VCF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1 (1-based)
            0x02, 0x00, 0x00, 0x00, // col_beg = 2 (1-based)
            0x00, 0x00, 0x00, 0x00, // col_end = None (1-based)
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x08, 0x00, 0x00, 0x00, // l_nm = 8
            b's', b'q', b'0', 0x00, // names[0] = "sq0"
            b's', b'q', b'1', 0x00, // names[1] = "sq1"
        ];
        let mut reader = &src[..];

        let actual = read_aux(&mut reader)?;

        let names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();
        let expected = crate::binning_index::index::header::Builder::vcf()
            .set_reference_sequence_names(names)
            .build();

        assert_eq!(actual, Some(expected));

        Ok(())
    }

    #[test]
    fn test_read_reference_sequence_name_index() -> Result<(), ReadError> {
        let data = [0x01, 0x00, 0x00, 0x00]; // col_seq = 1
        let mut reader = &data[..];
        assert_eq!(read_reference_sequence_name_index(&mut reader)?, 0);

        let data = [0xff, 0xff, 0xff, 0xff]; // col_seq = -1
        let mut reader = &data[..];
        assert!(matches!(
            read_reference_sequence_name_index(&mut reader),
            Err(ReadError::InvalidReferenceSequenceNameIndex(_))
        ));

        Ok(())
    }

    #[test]
    fn test_read_start_position_index() -> Result<(), ReadError> {
        let data = [0x04, 0x00, 0x00, 0x00]; // col_beg = 4
        let mut reader = &data[..];
        assert_eq!(read_start_position_index(&mut reader)?, 3);

        let data = [0xff, 0xff, 0xff, 0xff]; // col_beg = -1
        let mut reader = &data[..];
        assert!(matches!(
            read_start_position_index(&mut reader),
            Err(ReadError::InvalidStartPositionIndex(_))
        ));

        Ok(())
    }

    #[test]
    fn test_read_end_position_index() -> Result<(), ReadError> {
        let data = [0x00, 0x00, 0x00, 0x00]; // col_end = 0
        let mut reader = &data[..];
        assert!(read_end_position_index(&mut reader)?.is_none());

        let data = [0x05, 0x00, 0x00, 0x00]; // col_end = 5
        let mut reader = &data[..];
        assert_eq!(read_end_position_index(&mut reader)?, Some(4));

        let data = [0xff, 0xff, 0xff, 0xff]; // col_end = -1
        let mut reader = &data[..];
        assert!(matches!(
            read_end_position_index(&mut reader),
            Err(ReadError::InvalidEndPositionIndex(_))
        ));

        Ok(())
    }
}
