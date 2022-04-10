use crate::record::color::ParseError::Parse;
use std::fmt::{Display, Formatter};
use std::str::FromStr;
use std::{error, fmt, num};

/// A BED record color.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Color(u8, u8, u8);

impl fmt::Display for Color {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.0, self.1, self.2)
    }
}

/// An error returned when a raw BED record score fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    MissingR,
    MissingG,
    MissingB,
    /// The input failed to be parsed as an integer.
    Parse(num::ParseIntError),
}

impl error::Error for ParseError {}

impl Display for ParseError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => write!(f, "parse error: {}", e),
            Self::MissingR => f.write_str("missing r"),
            Self::MissingG => f.write_str("missing g"),
            Self::MissingB => f.write_str("missing b"),
        }
    }
}

impl FromStr for Color {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut args = s.split(",");
        let r = args
            .next()
            .ok_or(ParseError::MissingR)?
            .parse::<u8>()
            .map_err(ParseError::Parse)?;
        let g = args
            .next()
            .ok_or(ParseError::MissingG)?
            .parse::<u8>()
            .map_err(ParseError::Parse)?;
        let b = args
            .next()
            .ok_or(ParseError::MissingB)?
            .parse::<u8>()
            .map_err(ParseError::Parse)?;
        Result::Ok(Color(r, g, b))
    }
}
