mod reference_sequence_names;

use std::{
    error, fmt,
    io::{self, Read},
    num::{self, NonZero},
};

use self::reference_sequence_names::read_reference_sequence_names;
use crate::binning_index::index::{
    Header,
    header::{Format, format},
};

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
    /// The header start position index is invalid.
    InvalidStartPositionIndex(num::TryFromIntError),
    /// The header end position index is invalid.
    InvalidEndPositionIndex(num::TryFromIntError),
    /// The header end position index value is invalid.
    InvalidEndPositionIndexValue,
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
            Self::InvalidStartPositionIndex(_) => write!(f, "invalid start position index"),
            Self::InvalidEndPositionIndex(_) => write!(f, "invalid end position index"),
            Self::InvalidEndPositionIndexValue => write!(f, "invalid end position index value"),
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
    let l_aux =
        read_i32_le(reader).and_then(|n| u64::try_from(n).map_err(ReadError::InvalidAuxLength))?;

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
    let format =
        read_i32_le(reader).and_then(|n| Format::try_from(n).map_err(ReadError::InvalidFormat))?;

    let col_seq = read_reference_sequence_name_index(reader)?;
    let col_beg = read_start_position_index(reader)?;
    let col_end = read_end_position_index(reader, format, col_beg)?;
    let meta = read_line_comment_prefix(reader)?;

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
    crate::io::reader::num::read_i32_le(reader).map_err(ReadError::Io)
}

fn read_reference_sequence_name_index<R>(reader: &mut R) -> Result<usize, ReadError>
where
    R: Read,
{
    read_i32_le(reader).and_then(|i| {
        usize::try_from(i)
            .and_then(NonZero::try_from)
            .map(|n| n.get() - 1)
            .map_err(ReadError::InvalidReferenceSequenceNameIndex)
    })
}

fn read_start_position_index<R>(reader: &mut R) -> Result<usize, ReadError>
where
    R: Read,
{
    read_i32_le(reader).and_then(|i| {
        usize::try_from(i)
            .and_then(NonZero::try_from)
            .map(|n| n.get() - 1)
            .map_err(ReadError::InvalidStartPositionIndex)
    })
}

fn read_end_position_index<R>(
    reader: &mut R,
    format: Format,
    start_position_index: usize,
) -> Result<Option<usize>, ReadError>
where
    R: Read,
{
    const SPECIALIZED_END_VALUE: i32 = 0;

    let n = read_i32_le(reader)?;

    if matches!(format, Format::Sam | Format::Vcf) {
        if n == SPECIALIZED_END_VALUE {
            Ok(None)
        } else {
            Err(ReadError::InvalidEndPositionIndexValue)
        }
    } else {
        let i = usize::try_from(n)
            .and_then(NonZero::try_from)
            .map(|m| m.get() - 1)
            .map_err(ReadError::InvalidEndPositionIndex)?;

        if i == start_position_index {
            Ok(None)
        } else {
            Ok(Some(i))
        }
    }
}

fn read_line_comment_prefix<R>(reader: &mut R) -> Result<u8, ReadError>
where
    R: Read,
{
    read_i32_le(reader).and_then(|b| u8::try_from(b).map_err(ReadError::InvalidLineCommentPrefix))
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;
    use crate::binning_index::index::header::format::CoordinateSystem;

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
        let src = [0x01, 0x00, 0x00, 0x00];
        assert_eq!(read_reference_sequence_name_index(&mut &src[..])?, 0);

        let src = [0xff, 0xff, 0xff, 0xff];
        assert!(matches!(
            read_reference_sequence_name_index(&mut &src[..]),
            Err(ReadError::InvalidReferenceSequenceNameIndex(_))
        ));

        let src = [0x00, 0x00, 0x00, 0x00];
        assert!(matches!(
            read_reference_sequence_name_index(&mut &src[..]),
            Err(ReadError::InvalidReferenceSequenceNameIndex(_))
        ));

        Ok(())
    }

    #[test]
    fn test_read_start_position_index() -> Result<(), ReadError> {
        let src = [0x06, 0x00, 0x00, 0x00];
        assert_eq!(read_start_position_index(&mut &src[..])?, 5);

        let src = [0xff, 0xff, 0xff, 0xff];
        assert!(matches!(
            read_start_position_index(&mut &src[..]),
            Err(ReadError::InvalidStartPositionIndex(_))
        ));

        let src = [0x00, 0x00, 0x00, 0x00];
        assert!(matches!(
            read_start_position_index(&mut &src[..]),
            Err(ReadError::InvalidStartPositionIndex(_))
        ));

        Ok(())
    }

    #[test]
    fn test_read_end_position_index() -> Result<(), ReadError> {
        let src = [0x00, 0x00, 0x00, 0x00];
        assert!(read_end_position_index(&mut &src[..], Format::Sam, 5)?.is_none());

        let src = [0x00, 0x00, 0x00, 0x00];
        assert!(read_end_position_index(&mut &src[..], Format::Vcf, 5)?.is_none());

        let src = [0x09, 0x00, 0x00, 0x00];
        let format = Format::Generic(CoordinateSystem::Gff);
        assert_eq!(read_end_position_index(&mut &src[..], format, 5)?, Some(8));

        let src = [0x09, 0x00, 0x00, 0x00];
        let format = Format::Generic(CoordinateSystem::Gff);
        assert!(read_end_position_index(&mut &src[..], format, 8)?.is_none());

        let src = [0x09, 0x00, 0x00, 0x00];
        assert!(matches!(
            read_end_position_index(&mut &src[..], Format::Sam, 5),
            Err(ReadError::InvalidEndPositionIndexValue)
        ));

        let src = [0x09, 0x00, 0x00, 0x00];
        assert!(matches!(
            read_end_position_index(&mut &src[..], Format::Vcf, 5),
            Err(ReadError::InvalidEndPositionIndexValue)
        ));

        let src = [0xff, 0xff, 0xff, 0xff];
        let format = Format::Generic(CoordinateSystem::Gff);
        assert!(matches!(
            read_end_position_index(&mut &src[..], format, 5),
            Err(ReadError::InvalidEndPositionIndex(_))
        ));

        let src = [0x00, 0x00, 0x00, 0x00];
        let format = Format::Generic(CoordinateSystem::Gff);
        assert!(matches!(
            read_end_position_index(&mut &src[..], format, 5),
            Err(ReadError::InvalidEndPositionIndex(_))
        ));

        Ok(())
    }

    #[test]
    fn test_read_line_comment_prefix() -> Result<(), ReadError> {
        let src = [b'#', 0x00, 0x00, 0x00];
        assert_eq!(read_line_comment_prefix(&mut &src[..])?, b'#');

        let src = [0xff, 0xff, 0xff, 0xff];
        assert!(matches!(
            read_line_comment_prefix(&mut &src[..]),
            Err(ReadError::InvalidLineCommentPrefix(_))
        ));

        Ok(())
    }
}
