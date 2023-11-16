use std::{
    error, fmt,
    io::{self, Read},
    num, str,
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::index::{
    header::{format, ReferenceSequenceNames},
    Header,
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
    InvalidReferenceSequenceIndex(num::TryFromIntError),
    /// The header reference sequence index value is invalid.
    InvalidReferenceSequenceIndexValue,
    /// The header start position index is invalid.
    InvalidStartPositionIndex(num::TryFromIntError),
    /// The header start position index value is invalid.
    InvalidStartPositionIndexValue,
    /// The header end position index is invalid.
    InvalidEndPositionIndex(num::TryFromIntError),
    /// The header end position index value is invalid.
    InvalidEndPositionIndexValue,
    /// The header line comment prefix is invalid.
    InvalidLineCommentPrefix(num::TryFromIntError),
    /// The header line skip count is invalid.
    InvalidLineSkipCount(num::TryFromIntError),
    /// The header names length is invalid.
    InvalidNamesLength(num::TryFromIntError),
    /// A header name is duplicated.
    DuplicateName(String),
    /// The header names is invalid.
    InvalidNames,
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidAuxLength(e) => Some(e),
            Self::InvalidFormat(e) => Some(e),
            Self::InvalidReferenceSequenceIndex(e) => Some(e),
            Self::InvalidStartPositionIndex(e) => Some(e),
            Self::InvalidEndPositionIndex(e) => Some(e),
            Self::InvalidLineCommentPrefix(e) => Some(e),
            Self::InvalidLineSkipCount(e) => Some(e),
            Self::InvalidNamesLength(e) => Some(e),
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
            Self::InvalidReferenceSequenceIndex(_) => write!(f, "invalid reference sequence index"),
            Self::InvalidReferenceSequenceIndexValue => {
                write!(f, "invalid reference sequence index value")
            }
            Self::InvalidStartPositionIndex(_) => write!(f, "invalid start position index"),
            Self::InvalidStartPositionIndexValue => write!(f, "invalid start position index value"),
            Self::InvalidEndPositionIndex(_) => write!(f, "invalid end position index"),
            Self::InvalidEndPositionIndexValue => write!(f, "invalid end position index value"),
            Self::InvalidLineCommentPrefix(_) => write!(f, "invalid line comment prefix"),
            Self::InvalidLineSkipCount(_) => write!(f, "invalid line skip count"),
            Self::InvalidNamesLength(_) => write!(f, "invalid names length"),
            Self::DuplicateName(name) => write!(f, "duplicate name: {name}"),
            Self::InvalidNames => write!(f, "invalid names"),
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

pub(crate) fn read_header<R>(reader: &mut R) -> Result<Header, ReadError>
where
    R: Read,
{
    use crate::index::header::Format;

    let format =
        read_i32(reader).and_then(|n| Format::try_from(n).map_err(ReadError::InvalidFormat))?;

    let col_seq = read_i32(reader).and_then(|i| {
        usize::try_from(i)
            .map_err(ReadError::InvalidReferenceSequenceIndex)
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or(ReadError::InvalidReferenceSequenceIndexValue)
            })
    })?;

    let col_beg = read_i32(reader).and_then(|i| {
        usize::try_from(i)
            .map_err(ReadError::InvalidStartPositionIndex)
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or(ReadError::InvalidStartPositionIndexValue)
            })
    })?;

    let col_end = read_i32(reader).and_then(|i| match i {
        0 => Ok(None),
        _ => usize::try_from(i)
            .map_err(ReadError::InvalidEndPositionIndex)
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or(ReadError::InvalidEndPositionIndexValue)
            })
            .map(Some),
    })?;

    let meta = read_i32(reader)
        .and_then(|b| u8::try_from(b).map_err(ReadError::InvalidLineCommentPrefix))?;

    let skip =
        read_i32(reader).and_then(|n| u32::try_from(n).map_err(ReadError::InvalidLineSkipCount))?;

    let names = read_names(reader)?;

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

fn read_i32<R>(reader: &mut R) -> Result<i32, ReadError>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>().map_err(ReadError::Io)
}

fn read_names<R>(reader: &mut R) -> Result<ReferenceSequenceNames, ReadError>
where
    R: Read,
{
    let l_nm = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| usize::try_from(n).map_err(ReadError::InvalidNamesLength))?;

    let mut names = vec![0; l_nm];
    reader.read_exact(&mut names)?;

    parse_names(&names)
}

fn parse_names(mut src: &[u8]) -> Result<ReferenceSequenceNames, ReadError> {
    const NUL: u8 = 0x00;

    let mut names = ReferenceSequenceNames::new();

    while let Some(i) = src.iter().position(|&b| b == NUL) {
        let (raw_name, rest) = src.split_at(i);

        let name =
            str::from_utf8(raw_name).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if !names.insert(name.into()) {
            return Err(ReadError::DuplicateName(name.into()));
        }

        src = &rest[1..];
    }

    if src.is_empty() {
        Ok(names)
    } else {
        Err(ReadError::InvalidNames)
    }
}

#[cfg(test)]
mod tests {
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

        let names = [String::from("sq0"), String::from("sq1")]
            .into_iter()
            .collect();
        let expected = crate::index::header::Builder::vcf()
            .set_reference_sequence_names(names)
            .build();

        assert_eq!(actual, Some(expected));

        Ok(())
    }

    #[test]
    fn test_parse_names() -> Result<(), ReadError> {
        let data = b"sq0\x00sq1\x00";
        let actual = parse_names(&data[..])?;
        let expected: ReferenceSequenceNames = [String::from("sq0"), String::from("sq1")]
            .into_iter()
            .collect();
        assert_eq!(actual, expected);

        let data = b"";
        assert!(parse_names(&data[..])?.is_empty());

        Ok(())
    }

    #[test]
    fn test_parse_names_with_duplicate_name() {
        let data = b"sq0\x00sq0\x00";

        assert!(matches!(
            parse_names(data),
            Err(ReadError::DuplicateName(s)) if s == "sq0"
        ));
    }

    #[test]
    fn test_parse_names_with_trailing_data() {
        let data = b"sq0\x00sq1\x00sq2";
        assert!(matches!(parse_names(data), Err(ReadError::InvalidNames)));
    }
}
