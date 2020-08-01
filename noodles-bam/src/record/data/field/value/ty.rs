//! BAM record data field value type.

use std::{convert::TryFrom, error, fmt};

/// A BAM record data field value type.
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
        write!(f, "{}", char::from(*self))
    }
}

/// An error returned when a byte fails to convert to a BAM record data field value type.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid value")
    }
}

impl TryFrom<u8> for Type {
    type Error = TryFromByteError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            b'A' => Ok(Self::Char),
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
            _ => Err(TryFromByteError(b)),
        }
    }
}

impl From<Type> for char {
    fn from(ty: Type) -> char {
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
    fn test_try_from_byte_for_type() {
        assert_eq!(Type::try_from(b'A'), Ok(Type::Char));
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

        assert_eq!(Type::try_from(0), Err(TryFromByteError(0)));
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
