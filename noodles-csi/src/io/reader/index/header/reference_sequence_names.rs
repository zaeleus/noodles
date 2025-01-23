use std::{
    error, fmt,
    io::{self, BufRead, BufReader, Read},
    num,
};

use bstr::BString;
use byteorder::{LittleEndian, ReadBytesExt};

use crate::binning_index::index::header::ReferenceSequenceNames;

/// An error returned when CSI header reference sequence names fail to be read.
#[derive(Debug)]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
    /// A reference sequence name is duplicated.
    DuplicateName(BString),
    /// Expected EOF.
    ExpectedEof,
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
            Self::DuplicateName(name) => write!(f, "duplicate name: {name}"),
            Self::ExpectedEof => write!(f, "expected EOF"),
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

pub(super) fn read_reference_sequence_names<R>(
    reader: &mut R,
) -> Result<ReferenceSequenceNames, ReadError>
where
    R: Read,
{
    let l_nm = reader
        .read_i32::<LittleEndian>()
        .map_err(ReadError::Io)
        .and_then(|n| u64::try_from(n).map_err(ReadError::InvalidLength))?;

    let mut names_reader = BufReader::new(reader.take(l_nm));
    read_names(&mut names_reader)
}

fn read_names<R>(reader: &mut R) -> Result<ReferenceSequenceNames, ReadError>
where
    R: BufRead,
{
    let mut names = ReferenceSequenceNames::new();

    while let Some(name) = read_name(reader)? {
        if !names.insert(name.clone()) {
            return Err(ReadError::DuplicateName(name));
        }
    }

    Ok(names)
}

fn read_name<R>(reader: &mut R) -> Result<Option<BString>, ReadError>
where
    R: BufRead,
{
    const NUL: u8 = 0x00;

    let mut buf = BString::default();

    match reader.read_until(NUL, &mut buf)? {
        0 => Ok(None),
        _ => {
            if buf.ends_with(&[NUL]) {
                buf.pop();
                Ok(Some(buf))
            } else {
                Err(ReadError::ExpectedEof)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_reference_sequence_names() -> Result<(), ReadError> {
        let src = [
            0x08, 0x00, 0x00, 0x00, // l_nm = 8
            b's', b'q', b'0', 0x00, // names[0] = "sq0"
            b's', b'q', b'1', 0x00, // names[1] = "sq1"
        ];
        let mut reader = &src[..];

        let actual = read_reference_sequence_names(&mut reader)?;
        let expected: ReferenceSequenceNames = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();

        assert_eq!(actual, expected);

        let src = [0x00, 0x00, 0x00, 0x00]; // l_nm = 0
        let mut reader = &src[..];
        assert!(read_reference_sequence_names(&mut reader)?.is_empty());

        Ok(())
    }

    #[test]
    fn test_read_names_with_multiple_buffer_fills() -> Result<(), ReadError> {
        let src = b"sq0\x00sq1\x00";

        let mut reader = BufReader::with_capacity(2, &src[..]);

        let actual = read_names(&mut reader)?;
        let expected: ReferenceSequenceNames = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_names_with_duplicate_name() {
        let src = b"sq0\x00sq0\x00";
        let mut reader = &src[..];

        assert!(matches!(
            read_names(&mut reader),
            Err(ReadError::DuplicateName(s)) if s == "sq0"
        ));
    }

    #[test]
    fn test_read_names_with_trailing_data() {
        let src = b"sq0\x00sq1";
        let mut reader = &src[..];

        assert!(matches!(
            read_names(&mut reader),
            Err(ReadError::ExpectedEof)
        ));
    }
}
