use std::{error, fmt};

use crate::header::parser::record::value::map::format::Type;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid { actual: String },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid { actual } => write!(
                f,
                "invalid input: expected {{Integer, Float, Character, String}}, got {actual}"
            ),
        }
    }
}

pub(super) fn parse_type(s: &str) -> Result<Type, ParseError> {
    match s {
        "" => Err(ParseError::Empty),
        "Integer" => Ok(Type::Integer),
        "Float" => Ok(Type::Float),
        "Character" => Ok(Type::Character),
        "String" => Ok(Type::String),
        _ => Err(ParseError::Invalid { actual: s.into() }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_type() {
        assert_eq!(parse_type("Integer"), Ok(Type::Integer));
        assert_eq!(parse_type("Float"), Ok(Type::Float));
        assert_eq!(parse_type("Character"), Ok(Type::Character));
        assert_eq!(parse_type("String"), Ok(Type::String));

        assert_eq!(parse_type(""), Err(ParseError::Empty));
        assert_eq!(
            parse_type("ndls"),
            Err(ParseError::Invalid {
                actual: String::from("ndls")
            })
        );
    }
}
