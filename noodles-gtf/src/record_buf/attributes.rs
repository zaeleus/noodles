//! GTF record attributes.

pub mod entry;

pub use self::entry::Entry;

use std::{error, fmt, str::FromStr};

/// GTF record attributes.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(Vec<Entry>);

impl Attributes {
    /// Returns whether there are any entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record_buf::Attributes;
    /// let attributes = Attributes::default();
    /// assert!(attributes.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of entries.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record_buf::Attributes;
    /// let attributes = Attributes::default();
    /// assert_eq!(attributes.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns the value at the given key.
    ///
    /// This returns the first match of the key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf::record_buf::{attributes::Entry, Attributes};
    /// let attributes = Attributes::from(vec![Entry::new("id", "g0")]);
    /// assert_eq!(attributes.get("id"), Some("g0"));
    /// assert!(attributes.get("source").is_none());
    /// ```
    pub fn get<'a>(&'a self, key: &str) -> Option<&'a str> {
        self.0
            .iter()
            .find(|entry| entry.key() == key)
            .map(|entry| entry.value())
    }
}

impl AsRef<[Entry]> for Attributes {
    fn as_ref(&self) -> &[Entry] {
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
    /// The input has an invalid entry.
    InvalidEntry(entry::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidEntry(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidEntry(_) => write!(f, "invalid entry"),
        }
    }
}

impl FromStr for Attributes {
    type Err = ParseError;

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        use self::entry::parse_entry;

        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut entries = Vec::new();

        while !s.is_empty() {
            let entry = parse_entry(&mut s).map_err(ParseError::InvalidEntry)?;
            entries.push(entry);
        }

        Ok(Self(entries))
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
            r#"gene_id "g0""#.parse::<Attributes>(),
            Ok(Attributes::from(vec![Entry::new("gene_id", "g0")]))
        );

        assert_eq!(
            r#"gene_id "g0"; transcript_id "t0";"#.parse(),
            Ok(Attributes::from(vec![
                Entry::new("gene_id", "g0"),
                Entry::new("transcript_id", "t0")
            ]))
        );

        assert_eq!(
            r#"gene_id "g0";transcript_id "t0";"#.parse::<Attributes>(),
            Ok(Attributes::from(vec![
                Entry::new("gene_id", "g0"),
                Entry::new("transcript_id", "t0")
            ]))
        );

        assert_eq!(
            r#"gene_id "g0";  transcript_id "t0";"#.parse::<Attributes>(),
            Ok(Attributes::from(vec![
                Entry::new("gene_id", "g0"),
                Entry::new("transcript_id", "t0")
            ]))
        );

        assert_eq!("".parse::<Attributes>(), Err(ParseError::Empty));
        assert!(matches!(
            r#";"#.parse::<Attributes>(),
            Err(ParseError::InvalidEntry(_))
        ));
    }
}
