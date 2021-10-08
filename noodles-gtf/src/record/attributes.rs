//! GTF record attributes.

pub mod entry;

pub use self::entry::Entry;

use std::{error, fmt, ops::Deref, str::FromStr};

/// GTF record attributes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Attributes(Vec<Entry>);

impl Deref for Attributes {
    type Target = [Entry];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Entry>> for Attributes {
    fn from(entries: Vec<Entry>) -> Self {
        Self(entries)
    }
}

/// An error returned when raw GTF attributes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The input attributes has an invalid entry.
    InvalidEntry(entry::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidEntry(e) => write!(f, "invalid entry: {}", e),
        }
    }
}

impl FromStr for Attributes {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        const DELIMITER: &str = "; ";

        if s.is_empty() {
            Err(ParseError::Empty)
        } else if let Some(s) = s.strip_suffix(';') {
            s.split(DELIMITER)
                .map(|t| t.parse())
                .collect::<Result<Vec<_>, _>>()
                .map(Self::from)
                .map_err(ParseError::InvalidEntry)
        } else {
            Err(ParseError::Invalid)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(
            r#"gene_id "g0";"#.parse(),
            Ok(Attributes::from(vec![Entry::new("gene_id", "g0")]))
        );

        assert_eq!(
            r#"gene_id "g0"; transcript_id "t0";"#.parse(),
            Ok(Attributes::from(vec![
                Entry::new("gene_id", "g0"),
                Entry::new("transcript_id", "t0")
            ]))
        );

        assert_eq!("".parse::<Attributes>(), Err(ParseError::Empty));
        assert_eq!(
            r#"gene_id "g0""#.parse::<Attributes>(),
            Err(ParseError::Invalid)
        );
    }
}
