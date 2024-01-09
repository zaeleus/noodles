//! SAM record data field value type.

use std::fmt::{self, Write};

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
