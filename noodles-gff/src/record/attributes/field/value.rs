//! GFF record attributes field value.

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
    str::{self, FromStr},
};

const DELIMITER: char = ',';

/// A GFF record attribute field value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Value(Vec<String>);

impl Deref for Value {
    type Target = Vec<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Value {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use super::percent_encode;

        for (i, value) in self.iter().enumerate() {
            if i > 0 {
                DELIMITER.fmt(f)?;
            }

            percent_encode(value).fmt(f)?;
        }

        Ok(())
    }
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self(vec![String::from(s)])
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self(vec![s])
    }
}

impl From<Vec<String>> for Value {
    fn from(values: Vec<String>) -> Self {
        Self(values)
    }
}

/// An error returned when a raw GFF record attribute field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(str::Utf8Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Value {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use super::percent_decode;

        s.split(DELIMITER)
            .map(|s| {
                percent_decode(s)
                    .map(|t| t.into_owned())
                    .map_err(ParseError::Invalid)
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
        let value = Value::from("gene0");
        assert_eq!(value.to_string(), "gene0");

        let value = Value::from("13,21");
        assert_eq!(value.to_string(), "13%2C21");

        let value = Value::from(vec![String::from("13"), String::from("21")]);
        assert_eq!(value.to_string(), "13,21");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("gene0".parse(), Ok(Value::from("gene0")));
        assert_eq!("13%2C21".parse(), Ok(Value::from("13,21")));
        assert_eq!(
            "13,21".parse(),
            Ok(Value::from(vec![String::from("13"), String::from("21")]))
        );
    }
}
