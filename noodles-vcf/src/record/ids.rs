//! VCF record IDs.

use std::{error, fmt, ops::Deref, ops::DerefMut, str::FromStr};

use indexmap::IndexSet;

const DELIMITER: char = ';';

/// VCF record IDs (`ID`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Ids(IndexSet<String>);

impl Deref for Ids {
    type Target = IndexSet<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Ids {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Ids {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, id) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "{DELIMITER}")?;
            }

            f.write_str(id)?;
        }

        Ok(())
    }
}

/// An error returned when a raw VCF record ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The list of IDs has a duplicate.
    DuplicateId(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::DuplicateId(id) => write!(f, "duplicate ID: {id}"),
        }
    }
}

impl FromStr for Ids {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut ids = IndexSet::new();

        for raw_id in s.split(DELIMITER) {
            if !ids.insert(raw_id.into()) {
                return Err(ParseError::DuplicateId(raw_id.into()));
            }
        }

        Ok(Self(ids))
    }
}

impl Extend<String> for Ids {
    fn extend<T: IntoIterator<Item = String>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<String> for Ids {
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        let mut ids = Self::default();
        ids.extend(iter);
        ids
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert!(Ids::default().to_string().is_empty());

        let id0 = String::from("nd0");
        let id1 = String::from("nd1");

        let ids: Ids = [id0.clone()].into_iter().collect();
        assert_eq!(ids.to_string(), "nd0");

        let ids: Ids = [id0, id1].into_iter().collect();
        assert_eq!(ids.to_string(), "nd0;nd1");
    }

    #[test]
    fn test_from_str() {
        let id0 = String::from("nd0");
        let id1 = String::from("nd1");

        let expected: Ids = [id0.clone()].into_iter().collect();
        assert_eq!("nd0".parse(), Ok(expected));

        let expected: Ids = [id0, id1].into_iter().collect();
        assert_eq!("nd0;nd1".parse(), Ok(expected));

        assert_eq!("".parse::<Ids>(), Err(ParseError::Empty));
        assert_eq!(
            "nd0;nd0".parse::<Ids>(),
            Err(ParseError::DuplicateId(String::from("nd0")))
        );
    }
}
