//! VCF record IDs.

use std::{collections::HashSet, error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

const DELIMITER: char = ';';

/// VCF record IDs (`ID`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Ids(Vec<String>);

impl Deref for Ids {
    type Target = [String];

    fn deref(&self) -> &Self::Target {
        &self.0
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
    InvalidId(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::DuplicateId(id) => write!(f, "duplicate ID: {}", id),
            Self::InvalidId(s) => write!(f, "invalid ID: {}", s),
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
                let mut set: HashSet<String> = HashSet::new();
                let mut ids = Vec::new();

                for id in s.split(DELIMITER) {
                    if !set.insert(id.into()) {
                        return Err(ParseError::DuplicateId(id.into()));
                    } else if !is_valid_id(id) {
                        return Err(ParseError::InvalidId(id.into()));
                    }

                    ids.push(id.into());
                }

                Ok(Self(ids))
            }
        }
    }
}

fn is_valid_id(s: &str) -> bool {
    s.chars().all(|c| !c.is_ascii_whitespace())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Ids::default().to_string(), ".");
        assert_eq!(Ids(vec![String::from("nd0")]).to_string(), "nd0");
        assert_eq!(
            Ids(vec![String::from("nd0"), String::from("nd1")]).to_string(),
            "nd0;nd1"
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!(".".parse(), Ok(Ids::default()));
        assert_eq!("nd0".parse(), Ok(Ids(vec![String::from("nd0")])));
        assert_eq!(
            "nd0;nd1".parse(),
            Ok(Ids(vec![String::from("nd0"), String::from("nd1")]))
        );

        assert_eq!("".parse::<Ids>(), Err(ParseError::Empty));
        assert_eq!(
            "nd0;nd0".parse::<Ids>(),
            Err(ParseError::DuplicateId(String::from("nd0")))
        );
        assert_eq!(
            "nd 0".parse::<Ids>(),
            Err(ParseError::InvalidId(String::from("nd 0")))
        )
    }
}
