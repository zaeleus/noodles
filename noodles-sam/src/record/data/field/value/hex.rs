//! SAM record data field hex.

use std::{error, fmt, str::FromStr};

/// A SAM record data field hex value.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Hex(String);

impl AsRef<str> for Hex {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Hex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

/// An error returned when a raw hex string fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The length is invalid.
    InvalidLength {
        /// The actual length.
        actual: usize,
    },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidLength { actual } => {
                write!(
                    f,
                    "invalid length: expected length to be even, got {actual}"
                )
            }
        }
    }
}

impl FromStr for Hex {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from(s.as_bytes())
    }
}

impl TryFrom<&[u8]> for Hex {
    type Error = ParseError;

    fn try_from(buf: &[u8]) -> Result<Self, Self::Error> {
        if !is_even_length(buf.len()) {
            Err(ParseError::InvalidLength { actual: buf.len() })
        } else if buf.iter().copied().all(is_upper_ascii_hexdigit) {
            // SAFETY: `buf` is guaranteed to contain only hex digits.
            Ok(Self(String::from_utf8(buf.into()).unwrap()))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

fn is_even_length(n: usize) -> bool {
    n % 2 == 0
}

fn is_upper_ascii_hexdigit(n: u8) -> bool {
    matches!(n, b'0'..=b'9' | b'A'..=b'F')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_slice_for_hex() {
        assert_eq!(Hex::try_from(&b"CAFE"[..]), Ok(Hex(String::from("CAFE"))));

        assert_eq!(Hex::try_from(&b"cafe"[..]), Err(ParseError::Invalid));
        assert_eq!(Hex::try_from(&b"NDLS"[..]), Err(ParseError::Invalid));
        assert_eq!(
            Hex::try_from(&b"CAF"[..]),
            Err(ParseError::InvalidLength { actual: 3 })
        );
        assert_eq!(
            Hex::try_from(&b"CAFE0"[..]),
            Err(ParseError::InvalidLength { actual: 5 })
        );
    }
}
