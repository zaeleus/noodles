use std::{
    error, fmt,
    io::{self, Read},
    num, str,
};

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
    DuplicateName(String),
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
        .and_then(|n| usize::try_from(n).map_err(ReadError::InvalidLength))?;

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
        Err(ReadError::ExpectedEof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert!(matches!(parse_names(data), Err(ReadError::ExpectedEof)));
    }
}
