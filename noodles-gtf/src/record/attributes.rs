//! GTF record attributes.

pub mod entry;

pub use self::entry::Entry;

use std::{
    error,
    fmt::{self, Write},
    ops::Deref,
    str::FromStr,
};

const DELIMITER: char = ' ';

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

impl fmt::Display for Attributes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, entry) in self.0.iter().enumerate() {
            write!(f, "{}", entry)?;

            if i < self.0.len() - 1 {
                f.write_char(DELIMITER)?;
            }
        }

        Ok(())
    }
}

/// An error returned when raw GTF attributes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The input has an invalid entry.
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

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut entries = Vec::new();

        while let Some(i) = s.chars().position(|c| c == entry::TERMINATOR) {
            let (raw_entry, tail) = s.split_at(i + 1);

            let entry = raw_entry.parse().map_err(ParseError::InvalidEntry)?;
            entries.push(entry);

            s = consume_space(tail)?;
        }

        if s.is_empty() {
            Ok(Self(entries))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

// _GTF2.2: A Gene Annotation Format_ (2013-02-25): "Attributes must end in a semicolon which must
// then be separated from the start of any subsequent attribute by exactly one space character (NOT
// a tab character)."
fn consume_space(s: &str) -> Result<&str, ParseError> {
    if s.is_empty() {
        return Ok(s);
    } else if let Some(t) = s.strip_prefix(DELIMITER) {
        if t.strip_prefix(DELIMITER).is_none() {
            return Ok(t);
        }
    }

    Err(ParseError::Invalid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let attributes = Attributes::from(vec![Entry::new("gene_id", "g0")]);
        assert_eq!(attributes.to_string(), r#"gene_id "g0";"#);

        let attributes = Attributes::from(vec![
            Entry::new("gene_id", "g0"),
            Entry::new("transcript_id", "t0"),
        ]);
        assert_eq!(
            attributes.to_string(),
            r#"gene_id "g0"; transcript_id "t0";"#
        );
    }

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
        assert_eq!(
            r#"gene_id "g0";transcript_id "t0";"#.parse::<Attributes>(),
            Err(ParseError::Invalid)
        );
        assert_eq!(
            r#"gene_id "g0";  transcript_id "t0";"#.parse::<Attributes>(),
            Err(ParseError::Invalid)
        );
        assert!(matches!(
            r#";"#.parse::<Attributes>(),
            Err(ParseError::InvalidEntry(_))
        ));
    }
}
