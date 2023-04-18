use std::{error, fmt};

use crate::record::ReadName;

const MAX_LENGTH: usize = 254;

/// An error when a raw SAM record read name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is too long.
    ExceedsMaxLength(usize),
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::ExceedsMaxLength(actual) => {
                write!(f, "expected input to be < {MAX_LENGTH}, got {actual}")
            }
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub(crate) fn parse_read_name(src: &[u8]) -> Result<ReadName, ParseError> {
    if src.is_empty() {
        Err(ParseError::Empty)
    } else if src.len() > MAX_LENGTH {
        Err(ParseError::ExceedsMaxLength(src.len()))
    } else if is_valid_name(src) {
        Ok(ReadName(src.to_vec()))
    } else {
        Err(ParseError::Invalid)
    }
}

fn is_valid_name_char(b: u8) -> bool {
    matches!(b, b'!'..=b'?' | b'A'..=b'~')
}

fn is_valid_name(s: &[u8]) -> bool {
    s.iter().copied().all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_name() -> Result<(), crate::record::read_name::ParseError> {
        assert_eq!(parse_read_name(b"r0"), Ok(ReadName::try_new("r0")?));

        assert_eq!(parse_read_name(b""), Err(ParseError::Empty));
        assert_eq!(parse_read_name(b"r 0"), Err(ParseError::Invalid));
        assert_eq!(parse_read_name(b"@r0"), Err(ParseError::Invalid));

        let src = vec![b'n'; MAX_LENGTH + 1];
        assert_eq!(
            parse_read_name(&src),
            Err(ParseError::ExceedsMaxLength(src.len()))
        );

        Ok(())
    }
}
