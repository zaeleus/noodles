//! BAM record read name.

use std::{error, fmt, ops::Deref};

pub(crate) const MISSING: &[u8] = b"*";

// § 1.4 "The alignment section: mandatory fields" (2021-06-03)
const MAX_LENGTH: usize = 254;

/// A BAM record read name.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ReadName(Vec<u8>);

impl Deref for ReadName {
    type Target = Vec<u8>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// An error returned when a raw BAM record read name fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromBytesError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for TryFromBytesError {}

impl fmt::Display for TryFromBytesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl TryFrom<Vec<u8>> for ReadName {
    type Error = TryFromBytesError;

    fn try_from(buf: Vec<u8>) -> Result<Self, Self::Error> {
        if is_valid_name(&buf) {
            Ok(Self(buf))
        } else {
            Err(TryFromBytesError::Invalid)
        }
    }
}

impl From<ReadName> for Vec<u8> {
    fn from(read_name: ReadName) -> Self {
        read_name.0
    }
}

fn is_valid_name_char(b: u8) -> bool {
    b.is_ascii_graphic() && b != b'@'
}

fn is_valid_name(buf: &[u8]) -> bool {
    !buf.is_empty()
        && buf != MISSING
        && buf.len() <= MAX_LENGTH
        && buf.iter().copied().all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_vec_u8_for_read_name() {
        assert_eq!(
            ReadName::try_from(b"r0".to_vec()),
            Ok(ReadName(b"r0".to_vec()))
        );

        assert_eq!(
            ReadName::try_from(Vec::new()),
            Err(TryFromBytesError::Invalid)
        );
        assert_eq!(
            ReadName::try_from(b"*".to_vec()),
            Err(TryFromBytesError::Invalid)
        );
        assert_eq!(
            ReadName::try_from(b"r 0".to_vec()),
            Err(TryFromBytesError::Invalid)
        );
        assert_eq!(
            ReadName::try_from(b"@r0".to_vec()),
            Err(TryFromBytesError::Invalid)
        );

        let buf = b"n".repeat(MAX_LENGTH + 1);
        assert_eq!(ReadName::try_from(buf), Err(TryFromBytesError::Invalid));
    }
}
