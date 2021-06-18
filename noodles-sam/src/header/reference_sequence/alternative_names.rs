//! SAM header reference sequence alternative names.

use std::{error, fmt, ops::Deref, str::FromStr};

const DELIMITER: char = ',';

/// SAM header reference sequence alternative names (`AN`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternativeNames(Vec<String>);

impl Deref for AlternativeNames {
    type Target = [String];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for AlternativeNames {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, name) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "{}", DELIMITER)?;
            }

            write!(f, "{}", name)?;
        }

        Ok(())
    }
}

/// An error returned when raw SAM header reference sequence alternative names fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The raw names has an invalid name.
    InvalidName(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidName(s) => write!(f, "invalid name: {}", s),
        }
    }
}

impl FromStr for AlternativeNames {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use crate::record::reference_sequence_name::is_valid_name;

        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        s.split(DELIMITER)
            .map(|name| {
                if is_valid_name(name) {
                    Ok(name.into())
                } else {
                    Err(ParseError::InvalidName(name.into()))
                }
            })
            .collect::<Result<_, _>>()
            .map(Self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let alternative_names = AlternativeNames(vec![String::from("0")]);
        assert_eq!(alternative_names.to_string(), "0");

        let alternative_names = AlternativeNames(vec![String::from("0"), String::from("SQ.0")]);
        assert_eq!(alternative_names.to_string(), "0,SQ.0");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("0".parse(), Ok(AlternativeNames(vec![String::from("0")])));
        assert_eq!(
            "0,SQ.0".parse(),
            Ok(AlternativeNames(vec![
                String::from("0"),
                String::from("SQ.0")
            ]))
        );

        assert_eq!("".parse::<AlternativeNames>(), Err(ParseError::Empty));
        assert_eq!(
            ",".parse::<AlternativeNames>(),
            Err(ParseError::InvalidName(String::from("")))
        );
        assert_eq!(
            "0,*".parse::<AlternativeNames>(),
            Err(ParseError::InvalidName(String::from("*")))
        );
    }
}
