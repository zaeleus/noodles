//! GTF record attribute entry.

use std::{error, fmt, str::FromStr};

const SEPARATOR: char = ' ';
const DOUBLE_QUOTES: char = '"';
pub(super) const DELIMITER: char = ';';

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
    /// use noodles_gtf::record_buf::attributes::Entry;
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
    /// use noodles_gtf::record_buf::attributes::Entry;
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
    /// use noodles_gtf::record_buf::attributes::Entry;
    /// let entry = Entry::new("gene_id", "gene0");
    /// assert_eq!(entry.value(), "gene0");
    /// ```
    pub fn value(&self) -> &str {
        &self.value
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

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else {
            parse_entry(&mut s)
        }
    }
}

pub(super) fn parse_entry(s: &mut &str) -> Result<Entry, ParseError> {
    let key = parse_key(s)?;
    let value = parse_value(s)?;
    discard_delimiter(s);
    Ok(Entry::new(key, value))
}

fn parse_key<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    let Some(i) = s.find(SEPARATOR) else {
        return Err(ParseError::Invalid);
    };

    let (key, rest) = s.split_at(i);
    *s = &rest[1..];

    Ok(key)
}

fn parse_value<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    if let Some(rest) = s.strip_prefix(DOUBLE_QUOTES) {
        *s = rest;
        parse_string(s)
    } else {
        parse_raw_value(s)
    }
}

fn parse_string<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    if let Some(i) = s.find(DOUBLE_QUOTES) {
        let (t, rest) = s.split_at(i);
        *s = &rest[1..];
        Ok(t)
    } else {
        Err(ParseError::Invalid)
    }
}

fn parse_raw_value<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    if let Some(i) = s.find(DELIMITER) {
        let (t, rest) = s.split_at(i);
        *s = rest;
        Ok(t)
    } else {
        Ok(s)
    }
}

fn discard_delimiter(s: &mut &str) {
    *s = s.trim_start();

    if let Some(rest) = s.strip_prefix(DELIMITER) {
        *s = rest;
    }

    *s = s.trim_start();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(
            r#"gene_id "g0""#.parse::<Entry>(),
            Ok(Entry::new("gene_id", "g0"))
        );
        assert_eq!(
            r#"gene_ids "g0;g1""#.parse::<Entry>(),
            Ok(Entry::new("gene_ids", "g0;g1"))
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
