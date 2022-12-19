//! VCF record IDs.

pub mod id;

pub use self::id::Id;

use std::{error, fmt, ops::Deref, ops::DerefMut, str::FromStr};

use indexmap::IndexSet;

use super::MISSING_FIELD;

const DELIMITER: char = ';';

/// VCF record IDs (`ID`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Ids(IndexSet<Id>);

impl Deref for Ids {
    type Target = IndexSet<Id>;

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
        if self.is_empty() {
            write!(f, "{}", MISSING_FIELD)
        } else {
            for (i, id) in self.iter().enumerate() {
                if i > 0 {
                    write!(f, "{}", DELIMITER)?;
                }

                f.write_str(id)?;
            }

            Ok(())
        }
    }
}

/// An error returned when a raw VCF record ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The list of IDs has a duplicate.
    DuplicateId(String),
    /// An ID is invalid.
    InvalidId(id::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidId(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::DuplicateId(id) => write!(f, "duplicate ID: {}", id),
            Self::InvalidId(_) => f.write_str("invalid ID"),
        }
    }
}

impl FromStr for Ids {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Self::default()),
            _ => {
                let mut ids = IndexSet::new();

                for raw_id in s.split(DELIMITER) {
                    let id: Id = raw_id.parse().map_err(ParseError::InvalidId)?;

                    if !ids.insert(id) {
                        return Err(ParseError::DuplicateId(raw_id.into()));
                    }
                }

                Ok(Self(ids))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), id::ParseError> {
        assert_eq!(Ids::default().to_string(), ".");

        let id0: Id = "nd0".parse()?;
        let id1: Id = "nd1".parse()?;
        assert_eq!(Ids([id0.clone()].into_iter().collect()).to_string(), "nd0");
        assert_eq!(Ids([id0, id1].into_iter().collect()).to_string(), "nd0;nd1");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), id::ParseError> {
        let id0: Id = "nd0".parse()?;
        let id1: Id = "nd1".parse()?;

        assert_eq!(".".parse(), Ok(Ids::default()));
        assert_eq!("nd0".parse(), Ok(Ids([id0.clone()].into_iter().collect())));
        assert_eq!("nd0;nd1".parse(), Ok(Ids([id0, id1].into_iter().collect())));

        assert_eq!("".parse::<Ids>(), Err(ParseError::Empty));
        assert_eq!(
            "nd0;nd0".parse::<Ids>(),
            Err(ParseError::DuplicateId(String::from("nd0")))
        );
        assert!(matches!(
            "nd 0".parse::<Ids>(),
            Err(ParseError::InvalidId(_))
        ));
        assert!(matches!(
            ";nd0".parse::<Ids>(),
            Err(ParseError::InvalidId(_))
        ));
        assert!(matches!(
            "nd0;;nd1".parse::<Ids>(),
            Err(ParseError::InvalidId(_))
        ));
        assert!(matches!(
            "nd0;".parse::<Ids>(),
            Err(ParseError::InvalidId(_))
        ));

        Ok(())
    }
}
