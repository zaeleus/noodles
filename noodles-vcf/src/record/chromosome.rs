//! VCF record chromosome.

mod parser;

use std::{error, fmt, str::FromStr};

use super::MISSING_FIELD;

/// A VCF record chromosome (`CHROM`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Chromosome {
    /// A reference sequence name.
    Name(String),
    /// A symbol.
    Symbol(String),
}

impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Name(name) => f.write_str(name),
            Self::Symbol(symbol) => write!(f, "<{}>", symbol),
        }
    }
}

/// An error returned when a raw VCF record chromosome fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is missing (`.`).
    Missing,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Missing => f.write_str("missing input (`.`)"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Chromosome {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Err(ParseError::Missing),
            _ => parser::parse(s)
                .map(|(_, value)| match value {
                    parser::Value::Name(t) => Self::Name(t.into()),
                    parser::Value::Symbol(t) => Self::Symbol(t.into()),
                })
                .map_err(|_| ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Chromosome::Name(String::from("sq0")).to_string(), "sq0");
        assert_eq!(Chromosome::Symbol(String::from("sq0")).to_string(), "<sq0>");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("sq0".parse(), Ok(Chromosome::Name(String::from("sq0"))));
        assert_eq!("<sq0>".parse(), Ok(Chromosome::Symbol(String::from("sq0"))));

        assert_eq!("".parse::<Chromosome>(), Err(ParseError::Empty));
        assert_eq!(".".parse::<Chromosome>(), Err(ParseError::Missing));
    }
}
