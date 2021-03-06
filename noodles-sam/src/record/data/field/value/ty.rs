//! SAM record data field value type.

use std::{error, fmt, str::FromStr};

/// A SAM record data field value type.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    /// Character (`A`).
    Char,
    /// Integer (`i`).
    Int,
    /// Single-precision floating-point (`f`).
    Float,
    /// String (`Z`).
    String,
    /// Hex string (`H`).
    Hex,
    /// Array (`B`).
    Array,
}

impl fmt::Display for Type {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

/// An error returned when a raw SAM record data field value type fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid data field type: expected {{A, i, f, Z, H, B}}, got {}",
            self.0
        )
    }
}

impl FromStr for Type {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Self::Char),
            "i" => Ok(Self::Int),
            "f" => Ok(Self::Float),
            "Z" => Ok(Self::String),
            "H" => Ok(Self::Hex),
            "B" => Ok(Self::Array),
            _ => Err(ParseError(s.into())),
        }
    }
}

impl From<Type> for char {
    fn from(ty: Type) -> Self {
        match ty {
            Type::Char => 'A',
            Type::Int => 'i',
            Type::Float => 'f',
            Type::String => 'Z',
            Type::Hex => 'H',
            Type::Array => 'B',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Type::Char.to_string(), "A");
        assert_eq!(Type::Int.to_string(), "i");
        assert_eq!(Type::Float.to_string(), "f");
        assert_eq!(Type::String.to_string(), "Z");
        assert_eq!(Type::Hex.to_string(), "H");
        assert_eq!(Type::Array.to_string(), "B");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("A".parse(), Ok(Type::Char));
        assert_eq!("i".parse(), Ok(Type::Int));
        assert_eq!("f".parse(), Ok(Type::Float));
        assert_eq!("Z".parse(), Ok(Type::String));
        assert_eq!("H".parse(), Ok(Type::Hex));
        assert_eq!("B".parse(), Ok(Type::Array));

        assert_eq!("".parse::<Type>(), Err(ParseError(String::from(""))));
        assert_eq!("n".parse::<Type>(), Err(ParseError(String::from("n"))));
        assert_eq!(
            "noodles".parse::<Type>(),
            Err(ParseError(String::from("noodles")))
        );
    }

    #[test]
    fn test_from_type_for_char() {
        assert_eq!(char::from(Type::Char), 'A');
        assert_eq!(char::from(Type::Int), 'i');
        assert_eq!(char::from(Type::Float), 'f');
        assert_eq!(char::from(Type::String), 'Z');
        assert_eq!(char::from(Type::Hex), 'H');
        assert_eq!(char::from(Type::Array), 'B');
    }
}
