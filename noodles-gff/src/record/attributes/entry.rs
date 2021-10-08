//! GFF record attribute entry.

use std::{
    borrow::Cow,
    error, fmt,
    str::{self, FromStr},
};

use percent_encoding::{percent_decode_str, utf8_percent_encode, AsciiSet, CONTROLS};

const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS
    .add(b'\t')
    .add(b'\n')
    .add(b'\r')
    .add(b'%')
    .add(b';')
    .add(b'=')
    .add(b'&')
    .add(b',');

const SEPARATOR: char = '=';

/// A GFF record attribute entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Entry {
    key: String,
    value: String,
}

impl Entry {
    /// Creates a GFF record attribute.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::record::attributes::Entry;
    /// let entry = Entry::new("gene_name", "gene0");
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

    /// Returns the key of the entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::record::attributes::Entry;
    /// let entry = Entry::new("gene_name", "gene0");
    /// assert_eq!(entry.key(), "gene_name");
    /// ```
    pub fn key(&self) -> &str {
        &self.key
    }

    /// Returns the value of the entry.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::record::attributes::Entry;
    /// let entry = Entry::new("gene_name", "gene0");
    /// assert_eq!(entry.value(), "gene0");
    /// ```
    pub fn value(&self) -> &str {
        &self.value
    }
}

impl fmt::Display for Entry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}{}",
            percent_encode(self.key()),
            SEPARATOR,
            percent_encode(self.value())
        )
    }
}

/// An error returned when a raw GFF record attribute entry fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// The entry key is missing.
    MissingKey,
    /// The entry key is invalid.
    InvalidKey(str::Utf8Error),
    /// The entry value is missing.
    MissingValue,
    /// The entry value is invalid.
    InvalidValue(str::Utf8Error),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::MissingKey => f.write_str("missing key"),
            Self::InvalidKey(e) => write!(f, "invalid key: {}", e),
            Self::MissingValue => f.write_str("missing value"),
            Self::InvalidValue(e) => write!(f, "invalid value: {}", e),
        }
    }
}

impl FromStr for Entry {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        match s.split_once(SEPARATOR) {
            Some((k, v)) => {
                let key = parse_key(k)?;

                let value = if v.is_empty() {
                    return Err(ParseError::MissingValue);
                } else {
                    percent_decode(v).map_err(ParseError::InvalidValue)?
                };

                Ok(Self::new(key, value))
            }
            None => Err(ParseError::Invalid),
        }
    }
}

fn parse_key(s: &str) -> Result<Cow<'_, str>, ParseError> {
    if s.is_empty() {
        Err(ParseError::MissingKey)
    } else {
        percent_decode(s).map_err(ParseError::InvalidKey)
    }
}

fn percent_decode(s: &str) -> Result<Cow<'_, str>, str::Utf8Error> {
    percent_decode_str(s).decode_utf8()
}

fn percent_encode(s: &str) -> Cow<'_, str> {
    utf8_percent_encode(s, PERCENT_ENCODE_SET).into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let entry = Entry::new("gene_name", "gene0");
        assert_eq!(entry.to_string(), "gene_name=gene0");

        let entry = Entry::new("%s", "13,21");
        assert_eq!(entry.to_string(), "%25s=13%2C21");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "gene_name=gene0".parse(),
            Ok(Entry::new("gene_name", "gene0"))
        );
        assert_eq!("%25s=13%2C21".parse(), Ok(Entry::new("%s", "13,21")));

        assert_eq!("".parse::<Entry>(), Err(ParseError::Empty));
        assert_eq!("gene_name".parse::<Entry>(), Err(ParseError::Invalid));
        assert_eq!("=gene0".parse::<Entry>(), Err(ParseError::MissingKey));
        assert_eq!("gene_name=".parse::<Entry>(), Err(ParseError::MissingValue));
    }
}
