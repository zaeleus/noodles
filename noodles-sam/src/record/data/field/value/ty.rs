//! SAM record data field value type.

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

/// A SAM record data field value type.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    /// Character (`A`).
    Char,
    /// 8-bit integer (`c`).
    Int8,
    /// 8-bit unsigned integer (`C`).
    UInt8,
    /// 16-bit integer (`s`).
    Int16,
    /// 16-bit unsigned integer (`S`).
    UInt16,
    /// 32-bit integer (`i`).
    Int32,
    /// 32-bit unsigned integer (`I`).
    UInt32,
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
        f.write_char(char::from(*self))
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
            "invalid data field type: expected {{A, c, C, s, S, i, I, f, Z, H, B}}, got {}",
            self.0
        )
    }
}

impl FromStr for Type {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Self::Char),
            "c" => Ok(Self::Int8),
            "C" => Ok(Self::UInt8),
            "s" => Ok(Self::Int16),
            "S" => Ok(Self::UInt16),
            "i" => Ok(Self::Int32),
            "I" => Ok(Self::UInt32),
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
            Type::Int8 => 'c',
            Type::UInt8 => 'C',
            Type::Int16 => 's',
            Type::UInt16 => 'S',
            Type::Int32 => 'i',
            Type::UInt32 => 'I',
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
        assert_eq!(Type::Int8.to_string(), "c");
        assert_eq!(Type::UInt8.to_string(), "C");
        assert_eq!(Type::Int16.to_string(), "s");
        assert_eq!(Type::UInt16.to_string(), "S");
        assert_eq!(Type::Int32.to_string(), "i");
        assert_eq!(Type::UInt32.to_string(), "I");
        assert_eq!(Type::Float.to_string(), "f");
        assert_eq!(Type::String.to_string(), "Z");
        assert_eq!(Type::Hex.to_string(), "H");
        assert_eq!(Type::Array.to_string(), "B");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("A".parse(), Ok(Type::Char));
        assert_eq!("c".parse(), Ok(Type::Int8));
        assert_eq!("C".parse(), Ok(Type::UInt8));
        assert_eq!("s".parse(), Ok(Type::Int16));
        assert_eq!("S".parse(), Ok(Type::UInt16));
        assert_eq!("i".parse(), Ok(Type::Int32));
        assert_eq!("I".parse(), Ok(Type::UInt32));
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
        assert_eq!(char::from(Type::Int8), 'c');
        assert_eq!(char::from(Type::UInt8), 'C');
        assert_eq!(char::from(Type::Int16), 's');
        assert_eq!(char::from(Type::UInt16), 'S');
        assert_eq!(char::from(Type::Int32), 'i');
        assert_eq!(char::from(Type::UInt32), 'I');
        assert_eq!(char::from(Type::Float), 'f');
        assert_eq!(char::from(Type::String), 'Z');
        assert_eq!(char::from(Type::Hex), 'H');
        assert_eq!(char::from(Type::Array), 'B');
    }
}
