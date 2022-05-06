//! SAM record data field value subtype.

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

/// A SAM record data field value subtype.
///
/// Only arrays have subtypes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Subtype {
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
}

impl fmt::Display for Subtype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char(char::from(*self))
    }
}

/// An error returned when a raw SAM record data field value subtype fails to parse.
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

impl FromStr for Subtype {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "c" => Ok(Self::Int8),
            "C" => Ok(Self::UInt8),
            "s" => Ok(Self::Int16),
            "S" => Ok(Self::UInt16),
            "i" => Ok(Self::Int32),
            "I" => Ok(Self::UInt32),
            "f" => Ok(Self::Float),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl TryFrom<u8> for Subtype {
    type Error = ParseError;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            b'c' => Ok(Self::Int8),
            b'C' => Ok(Self::UInt8),
            b's' => Ok(Self::Int16),
            b'S' => Ok(Self::UInt16),
            b'i' => Ok(Self::Int32),
            b'I' => Ok(Self::UInt32),
            b'f' => Ok(Self::Float),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl From<Subtype> for char {
    fn from(subtype: Subtype) -> Self {
        char::from(u8::from(subtype))
    }
}

impl From<Subtype> for u8 {
    fn from(subtype: Subtype) -> Self {
        match subtype {
            Subtype::Int8 => b'c',
            Subtype::UInt8 => b'C',
            Subtype::Int16 => b's',
            Subtype::UInt16 => b'S',
            Subtype::Int32 => b'i',
            Subtype::UInt32 => b'I',
            Subtype::Float => b'f',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Subtype::Int8.to_string(), "c");
        assert_eq!(Subtype::UInt8.to_string(), "C");
        assert_eq!(Subtype::Int16.to_string(), "s");
        assert_eq!(Subtype::UInt16.to_string(), "S");
        assert_eq!(Subtype::Int32.to_string(), "i");
        assert_eq!(Subtype::UInt32.to_string(), "I");
        assert_eq!(Subtype::Float.to_string(), "f");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("c".parse(), Ok(Subtype::Int8));
        assert_eq!("C".parse(), Ok(Subtype::UInt8));
        assert_eq!("s".parse(), Ok(Subtype::Int16));
        assert_eq!("S".parse(), Ok(Subtype::UInt16));
        assert_eq!("i".parse(), Ok(Subtype::Int32));
        assert_eq!("I".parse(), Ok(Subtype::UInt32));
        assert_eq!("f".parse(), Ok(Subtype::Float));

        assert_eq!("".parse::<Subtype>(), Err(ParseError::Invalid));
        assert_eq!("n".parse::<Subtype>(), Err(ParseError::Invalid));
        assert_eq!("noodles".parse::<Subtype>(), Err(ParseError::Invalid));
    }

    #[test]
    fn test_try_from_u8_for_subtype() {
        assert_eq!(Subtype::try_from(b'c'), Ok(Subtype::Int8));
        assert_eq!(Subtype::try_from(b'C'), Ok(Subtype::UInt8));
        assert_eq!(Subtype::try_from(b's'), Ok(Subtype::Int16));
        assert_eq!(Subtype::try_from(b'S'), Ok(Subtype::UInt16));
        assert_eq!(Subtype::try_from(b'i'), Ok(Subtype::Int32));
        assert_eq!(Subtype::try_from(b'I'), Ok(Subtype::UInt32));
        assert_eq!(Subtype::try_from(b'f'), Ok(Subtype::Float));

        assert_eq!(Subtype::try_from(b'n'), Err(ParseError::Invalid));
    }

    #[test]
    fn test_from_subtype_for_char() {
        assert_eq!(char::from(Subtype::Int8), 'c');
        assert_eq!(char::from(Subtype::UInt8), 'C');
        assert_eq!(char::from(Subtype::Int16), 's');
        assert_eq!(char::from(Subtype::UInt16), 'S');
        assert_eq!(char::from(Subtype::Int32), 'i');
        assert_eq!(char::from(Subtype::UInt32), 'I');
        assert_eq!(char::from(Subtype::Float), 'f');
    }

    #[test]
    fn test_from_subtype_for_u8() {
        assert_eq!(u8::from(Subtype::Int8), b'c');
        assert_eq!(u8::from(Subtype::UInt8), b'C');
        assert_eq!(u8::from(Subtype::Int16), b's');
        assert_eq!(u8::from(Subtype::UInt16), b'S');
        assert_eq!(u8::from(Subtype::Int32), b'i');
        assert_eq!(u8::from(Subtype::UInt32), b'I');
        assert_eq!(u8::from(Subtype::Float), b'f');
    }
}
