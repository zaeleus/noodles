//! SAM header reference sequence alternative locus.

use std::{error, fmt, str::FromStr};

const UNKNOWN: &str = "*";

/// A SAM header reference sequence alternative locus (`AH`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AlternativeLocus {
    /// A region in the primary assembly.
    Region(String, Option<(i32, i32)>),
    /// The region is unknown.
    Unknown,
}

impl fmt::Display for AlternativeLocus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Region(reference_sequence_name, Some((start, end))) => {
                write!(f, "{}:{}-{}", reference_sequence_name, start, end)
            }
            Self::Region(reference_sequence_name, None) => f.write_str(reference_sequence_name),
            Self::Unknown => f.write_str(UNKNOWN),
        }
    }
}

/// An error returned when raw SAM header reference sequence alternative locus fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The reference sequence name is invalid.
    InvalidReferenceSequenceName,
    /// The interval is invalid.
    InvalidInterval,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingReferenceSequenceName => f.write_str("missing reference sequence name"),
            Self::InvalidReferenceSequenceName => f.write_str("invalid reference sequence name"),
            Self::InvalidInterval => f.write_str("invalid interval"),
        }
    }
}

impl FromStr for AlternativeLocus {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use crate::record::reference_sequence_name::is_valid_name;

        match s {
            "" => Err(ParseError::Empty),
            UNKNOWN => Ok(Self::Unknown),
            _ => {
                let mut components = s.splitn(2, ':');

                let reference_sequence_name = components
                    .next()
                    .ok_or(ParseError::MissingReferenceSequenceName)?;

                if !is_valid_name(reference_sequence_name) {
                    return Err(ParseError::InvalidReferenceSequenceName);
                }

                let interval = components
                    .next()
                    .map(|t| parse_interval(t).map(Some))
                    .unwrap_or(Ok(None))?;

                Ok(Self::Region(reference_sequence_name.into(), interval))
            }
        }
    }
}

fn parse_interval(s: &str) -> Result<(i32, i32), ParseError> {
    let mut components = s.splitn(2, '-');

    let start = components
        .next()
        .ok_or(ParseError::InvalidInterval)
        .and_then(|t| t.parse().map_err(|_| ParseError::InvalidInterval))?;

    let end = components
        .next()
        .ok_or(ParseError::InvalidInterval)
        .and_then(|t| t.parse().map_err(|_| ParseError::InvalidInterval))?;

    Ok((start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(AlternativeLocus::Unknown.to_string(), "*");

        assert_eq!(
            AlternativeLocus::Region(String::from("sq0"), None).to_string(),
            "sq0"
        );

        assert_eq!(
            AlternativeLocus::Region(String::from("sq0"), Some((8, 13))).to_string(),
            "sq0:8-13"
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!("*".parse(), Ok(AlternativeLocus::Unknown));

        assert_eq!(
            "sq0".parse(),
            Ok(AlternativeLocus::Region(String::from("sq0"), None))
        );

        assert_eq!(
            "sq0:8-13".parse(),
            Ok(AlternativeLocus::Region(String::from("sq0"), Some((8, 13))))
        );

        assert_eq!("".parse::<AlternativeLocus>(), Err(ParseError::Empty));

        assert_eq!(
            "=".parse::<AlternativeLocus>(),
            Err(ParseError::InvalidReferenceSequenceName)
        );
    }
}
