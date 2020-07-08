use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

const DELIMITER: char = ';';

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
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
        }
    }
}

impl FromStr for Ids {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Self::default()),
            _ => Ok(Self(s.split(DELIMITER).map(|s| s.into()).collect())),
        }
    }
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
    }
}
