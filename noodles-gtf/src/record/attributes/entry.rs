//! GTF record attribute entry.

use std::{error, fmt, str::FromStr};

const SEPARATOR: char = ' ';
pub(super) const TERMINATOR: char = ';';

/// A GTF record attribute entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Entry {
    key: String,
    value: String,
}

impl Entry {
    /// Creates a GTF record attribute entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record::attributes::Entry;
    /// let entry = Entry::new("gene_id", "gene0");
    /// ```
    pub fn new<K, V>(key: K, value: V) -> Self
    where
        K: Into<String>,
        V: Into<String>,
    {
        Self {
            key: key.into(),
            value: value.into(),
        }
    }

    /// Creates a GTF record attribute entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record::attributes::Entry;
    /// let entry = Entry::new("gene_id", "gene0");
    /// assert_eq!(entry.key(), "gene_id");
    /// ```
    pub fn key(&self) -> &str {
        &self.key
    }

    /// Creates a GTF record attribute entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record::attributes::Entry;
    /// let entry = Entry::new("gene_id", "gene0");
    /// assert_eq!(entry.value(), "gene0");
    /// ```
    pub fn value(&self) -> &str {
        &self.value
    }
}

impl fmt::Display for Entry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, r#"{}{}"{}""#, self.key(), SEPARATOR, self.value())
    }
}

/// An error returned when a raw GTF record attribute entry fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl FromStr for Entry {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else {
            parse_entry(s)
        }
    }
}

fn parse_entry(s: &str) -> Result<Entry, ParseError> {
    match s.split_once(SEPARATOR) {
        Some((k, v)) => {
            let value = parse_value(v);
            Ok(Entry::new(k, value))
        }
        None => Err(ParseError::Invalid),
    }
}

fn parse_value(s: &str) -> String {
    s.trim_matches('"').into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let entry = Entry::new("gene_id", "g0");
        assert_eq!(entry.to_string(), r#"gene_id "g0""#);
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            r#"gene_id "g0""#.parse::<Entry>(),
            Ok(Entry::new("gene_id", "g0"))
        );
        assert_eq!(
            r#"gene_id """#.parse::<Entry>(),
            Ok(Entry::new("gene_id", ""))
        );
        assert_eq!(
            r#"gene_id 0"#.parse::<Entry>(),
            Ok(Entry::new("gene_id", "0"))
        );

        assert_eq!("".parse::<Entry>(), Err(ParseError::Empty));
        assert_eq!("gene_id".parse::<Entry>(), Err(ParseError::Invalid));
        assert_eq!(r#""""#.parse::<Entry>(), Err(ParseError::Invalid));
    }
}
