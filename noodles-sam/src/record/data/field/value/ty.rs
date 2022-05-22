//! SAM record data field value type.

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

/// A SAM record data field value type.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Type {
    /// Character (`A`).
    Character,
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
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Type {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Self::Character),
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
            _ => Err(ParseError::Invalid),
        }
    }
}

impl TryFrom<u8> for Type {
    type Error = ParseError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            b'A' => Ok(Self::Character),
            b'c' => Ok(Self::Int8),
            b'C' => Ok(Self::UInt8),
            b's' => Ok(Self::Int16),
            b'S' => Ok(Self::UInt16),
            b'i' => Ok(Self::Int32),
            b'I' => Ok(Self::UInt32),
            b'f' => Ok(Self::Float),
            b'Z' => Ok(Self::String),
            b'H' => Ok(Self::Hex),
            b'B' => Ok(Self::Array),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl From<Type> for char {
    fn from(ty: Type) -> Self {
        Self::from(u8::from(ty))
    }
}

impl From<Type> for u8 {
    fn from(ty: Type) -> Self {
        match ty {
            Type::Character => b'A',
            Type::Int8 => b'c',
            Type::UInt8 => b'C',
            Type::Int16 => b's',
            Type::UInt16 => b'S',
            Type::Int32 => b'i',
            Type::UInt32 => b'I',
            Type::Float => b'f',
            Type::String => b'Z',
            Type::Hex => b'H',
            Type::Array => b'B',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Type::Character.to_string(), "A");
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
        assert_eq!("A".parse(), Ok(Type::Character));
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

        assert_eq!("".parse::<Type>(), Err(ParseError::Invalid));
        assert_eq!("n".parse::<Type>(), Err(ParseError::Invalid));
        assert_eq!("noodles".parse::<Type>(), Err(ParseError::Invalid));
    }

    #[test]
    fn test_try_from_u8_for_type() {
        assert_eq!(Type::try_from(b'A'), Ok(Type::Character));
        assert_eq!(Type::try_from(b'c'), Ok(Type::Int8));
        assert_eq!(Type::try_from(b'C'), Ok(Type::UInt8));
        assert_eq!(Type::try_from(b's'), Ok(Type::Int16));
        assert_eq!(Type::try_from(b'S'), Ok(Type::UInt16));
        assert_eq!(Type::try_from(b'i'), Ok(Type::Int32));
        assert_eq!(Type::try_from(b'I'), Ok(Type::UInt32));
        assert_eq!(Type::try_from(b'f'), Ok(Type::Float));
        assert_eq!(Type::try_from(b'Z'), Ok(Type::String));
        assert_eq!(Type::try_from(b'H'), Ok(Type::Hex));
        assert_eq!(Type::try_from(b'B'), Ok(Type::Array));

        assert_eq!(Type::try_from(b'n'), Err(ParseError::Invalid));
    }

    #[test]
    fn test_from_type_for_char() {
        assert_eq!(char::from(Type::Character), 'A');
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

    #[test]
    fn test_from_type_for_u8() {
        assert_eq!(u8::from(Type::Character), b'A');
        assert_eq!(u8::from(Type::Int8), b'c');
        assert_eq!(u8::from(Type::UInt8), b'C');
        assert_eq!(u8::from(Type::Int16), b's');
        assert_eq!(u8::from(Type::UInt16), b'S');
        assert_eq!(u8::from(Type::Int32), b'i');
        assert_eq!(u8::from(Type::UInt32), b'I');
        assert_eq!(u8::from(Type::Float), b'f');
        assert_eq!(u8::from(Type::String), b'Z');
        assert_eq!(u8::from(Type::Hex), b'H');
        assert_eq!(u8::from(Type::Array), b'B');
    }
}
