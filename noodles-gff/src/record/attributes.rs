//! GFF record attributes and entry.

pub mod entry;

pub use self::entry::Entry;

use std::{ops::Deref, str::FromStr};

const DELIMITER: char = ';';

/// GFF record attributes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(Vec<Entry>);

impl Deref for Attributes {
    type Target = [Entry];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// An error returned when raw attributes fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input attributes has an invalid entry.
    InvalidEntry(entry::ParseError),
}

impl From<Vec<Entry>> for Attributes {
    fn from(entries: Vec<Entry>) -> Self {
        Self(entries)
    }
}

impl FromStr for Attributes {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        s.split(DELIMITER)
            .map(|t| t.parse())
            .collect::<Result<Vec<_>, _>>()
            .map(Self::from)
            .map_err(ParseError::InvalidEntry)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = "gene_id=ndls0;gene_name=gene0";
        let actual = s.parse::<Attributes>()?;
        let expected = Attributes::from(vec![
            Entry::new(String::from("gene_id"), String::from("ndls0")),
            Entry::new(String::from("gene_name"), String::from("gene0")),
        ]);
        assert_eq!(actual, expected);

        let s = "gene_id=ndls0";
        let actual = s.parse::<Attributes>()?;
        let expected = Attributes::from(vec![Entry::new(
            String::from("gene_id"),
            String::from("ndls0"),
        )]);
        assert_eq!(actual, expected);

        let actual = "".parse::<Attributes>()?;
        let expected = Attributes::default();
        assert_eq!(actual, expected);

        Ok(())
    }
}
