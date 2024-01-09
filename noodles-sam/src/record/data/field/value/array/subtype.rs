//! SAM record data field value subtype.

use std::fmt::{self, Write};

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
