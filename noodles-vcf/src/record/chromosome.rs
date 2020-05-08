mod parser;

use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Chromosome {
    Name(String),
    Reference(String),
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid chromosome: {}", self.0)
    }
}

impl FromStr for Chromosome {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError(s.into()));
        }

        let (_, value) = parser::parse(s).map_err(|_| ParseError(s.into()))?;

        match value {
            parser::Value::Name(s) => Ok(Self::Name(s.into())),
            parser::Value::Reference(s) => Ok(Self::Reference(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!(
            "sq0".parse::<Chromosome>()?,
            Chromosome::Name(String::from("sq0"))
        );

        assert_eq!(
            "<sq0>".parse::<Chromosome>()?,
            Chromosome::Reference(String::from("sq0"))
        );

        assert!("".parse::<Chromosome>().is_err());

        Ok(())
    }
}
