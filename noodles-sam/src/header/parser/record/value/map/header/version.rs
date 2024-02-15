use std::{error, fmt};

use crate::header::record::value::map::header::Version;

pub(crate) fn parse_version(src: &[u8]) -> Result<Version, ParseError> {
    const DELIMITER: u8 = b'.';

    fn split_once(buf: &[u8], delimiter: u8) -> Option<(&[u8], &[u8])> {
        let i = buf.iter().position(|&b| b == delimiter)?;
        Some((&buf[..i], &buf[i + 1..]))
    }

    match split_once(src, DELIMITER) {
        Some((a, b)) => {
            let major = lexical_core::parse(a).map_err(ParseError::InvalidMajorVersion)?;
            let minor = lexical_core::parse(b).map_err(ParseError::InvalidMinorVersion)?;
            Ok(Version::new(major, minor))
        }
        None => Err(ParseError::Invalid),
    }
}

/// An error returned when a SAM header header version fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The major version is invalid.
    InvalidMajorVersion(lexical_core::Error),
    /// The minor version is invalid.
    InvalidMinorVersion(lexical_core::Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidMajorVersion(e) | Self::InvalidMinorVersion(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidMajorVersion(_) => write!(f, "invalid major version"),
            Self::InvalidMinorVersion(_) => write!(f, "invalid minor version"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_version() {
        assert_eq!(parse_version(b"1.6"), Ok(Version::new(1, 6)));

        assert_eq!(parse_version(b""), Err(ParseError::Invalid));
        assert_eq!(parse_version(b"1"), Err(ParseError::Invalid));

        assert!(matches!(
            parse_version(b"."),
            Err(ParseError::InvalidMajorVersion(_))
        ));

        assert!(matches!(
            parse_version(b"x.6"),
            Err(ParseError::InvalidMajorVersion(_))
        ));

        assert!(matches!(
            parse_version(b"1.x"),
            Err(ParseError::InvalidMinorVersion(_))
        ));

        assert!(matches!(
            parse_version(b"1.6.1"),
            Err(ParseError::InvalidMinorVersion(_))
        ));
    }
}
