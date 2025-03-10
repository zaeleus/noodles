//! GTF record frame.

use std::{error, fmt, num, str::FromStr};

const MIN: u8 = 0;
const MAX: u8 = 2;

/// A GTF record frame.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Frame(u8);

impl fmt::Display for Frame {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// An error returned when a raw GTF record frame fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(num::ParseIntError),
    /// The value is invalid.
    InvalidValue(u8),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(_) => f.write_str("invalid input"),
            Self::InvalidValue(n) => {
                write!(f, "invalid value: expected {MIN}..={MAX}, got {n}")
            }
        }
    }
}

impl FromStr for Frame {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else {
            s.parse::<u8>()
                .map_err(ParseError::Invalid)
                .and_then(Self::try_from)
        }
    }
}

impl TryFrom<u8> for Frame {
    type Error = ParseError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        if n <= MAX {
            Ok(Self(n))
        } else {
            Err(ParseError::InvalidValue(n))
        }
    }
}

impl From<Frame> for u8 {
    fn from(frame: Frame) -> Self {
        frame.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Frame(0).to_string(), "0");
        assert_eq!(Frame(1).to_string(), "1");
        assert_eq!(Frame(2).to_string(), "2");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("0".parse(), Ok(Frame(0)));
        assert_eq!("1".parse(), Ok(Frame(1)));
        assert_eq!("2".parse(), Ok(Frame(2)));

        assert_eq!("".parse::<Frame>(), Err(ParseError::Empty));
        assert!(matches!("n".parse::<Frame>(), Err(ParseError::Invalid(_))));
        assert_eq!("3".parse::<Frame>(), Err(ParseError::InvalidValue(3)));
    }

    #[test]
    fn test_try_from_u8_for_frame() {
        assert_eq!(Frame::try_from(0), Ok(Frame(0)));
        assert_eq!(Frame::try_from(1), Ok(Frame(1)));
        assert_eq!(Frame::try_from(2), Ok(Frame(2)));
        assert_eq!(Frame::try_from(3), Err(ParseError::InvalidValue(3)));
    }

    #[test]
    fn test_from_frame_for_u8() {
        assert_eq!(u8::from(Frame(0)), 0);
    }
}
