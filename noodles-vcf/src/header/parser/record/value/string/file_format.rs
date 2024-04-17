use std::{error, fmt};

use crate::header::FileFormat;

/// An error returned when a VCF header file format value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The prefix is invalid.
    InvalidPrefix,
    /// The version is invalid.
    InvalidVersion,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPrefix => write!(f, "invalid prefix"),
            Self::InvalidVersion => write!(f, "invalid version"),
        }
    }
}

pub fn parse_file_format(mut src: &[u8]) -> Result<FileFormat, ParseError> {
    const PREFIX: &[u8] = b"VCFv";
    const DELIMITER: u8 = b'.';

    fn split_once(buf: &[u8], delimiter: u8) -> Option<(&[u8], &[u8])> {
        let i = buf.iter().position(|&b| b == delimiter)?;
        Some((&buf[..i], &buf[i + 1..]))
    }

    src = src.strip_prefix(PREFIX).ok_or(ParseError::InvalidPrefix)?;

    match split_once(src, DELIMITER) {
        Some((a, b)) => {
            let major = parse_u32(a)?;
            let minor = parse_u32(b)?;
            Ok(FileFormat::new(major, minor))
        }
        None => Err(ParseError::InvalidVersion),
    }
}

fn parse_u32(src: &[u8]) -> Result<u32, ParseError> {
    const OFFSET: u8 = b'0';

    src.iter().try_fold(0, |n: u32, b| {
        let d = match b.checked_sub(OFFSET) {
            Some(n @ 0..=9) => u32::from(n),
            _ => return Err(ParseError::InvalidVersion),
        };

        n.checked_mul(10)
            .and_then(|m| m.checked_add(d))
            .ok_or(ParseError::InvalidVersion)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_file_format() {
        assert_eq!(parse_file_format(b"VCFv4.3"), Ok(FileFormat::new(4, 3)));

        assert_eq!(parse_file_format(b""), Err(ParseError::InvalidPrefix));
        assert_eq!(parse_file_format(b"4.3"), Err(ParseError::InvalidPrefix));
        assert_eq!(
            parse_file_format(b"NDLv4.3"),
            Err(ParseError::InvalidPrefix)
        );

        assert_eq!(parse_file_format(b"VCFv"), Err(ParseError::InvalidVersion));
        assert_eq!(parse_file_format(b"VCFvx"), Err(ParseError::InvalidVersion));
        assert_eq!(parse_file_format(b"VCFv4"), Err(ParseError::InvalidVersion));
        assert_eq!(
            parse_file_format(b"VCFvx.3"),
            Err(ParseError::InvalidVersion)
        );
        assert_eq!(
            parse_file_format(b"VCFvx.4"),
            Err(ParseError::InvalidVersion)
        );
    }
}
