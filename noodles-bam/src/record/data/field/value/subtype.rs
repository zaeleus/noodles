//! BAM record data field value subtype.

use std::{error, fmt};

/// A BAM record data field value subtype.
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
        write!(f, "{}", char::from(*self))
    }
}

/// An error returned when a byte fails to convert to a BAM record data field value subtype.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteError(u8);

impl error::Error for TryFromByteError {}

impl fmt::Display for TryFromByteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expected {{0x63, 0x43, 0x73, 0x53, 0x69, 0x49, 0x66}}, got 0x{:02x}",
            self.0
        )
    }
}

impl TryFrom<u8> for Subtype {
    type Error = TryFromByteError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            b'c' => Ok(Self::Int8),
            b'C' => Ok(Self::UInt8),
            b's' => Ok(Self::Int16),
            b'S' => Ok(Self::UInt16),
            b'i' => Ok(Self::Int32),
            b'I' => Ok(Self::UInt32),
            b'f' => Ok(Self::Float),
            _ => Err(TryFromByteError(b)),
        }
    }
}

impl From<Subtype> for char {
    fn from(subtype: Subtype) -> Self {
        match subtype {
            Subtype::Int8 => 'c',
            Subtype::UInt8 => 'C',
            Subtype::Int16 => 's',
            Subtype::UInt16 => 'S',
            Subtype::Int32 => 'i',
            Subtype::UInt32 => 'I',
            Subtype::Float => 'f',
        }
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
    fn test_try_from_char_for_subtype() {
        assert_eq!(Subtype::try_from(b'c'), Ok(Subtype::Int8));
        assert_eq!(Subtype::try_from(b'C'), Ok(Subtype::UInt8));
        assert_eq!(Subtype::try_from(b's'), Ok(Subtype::Int16));
        assert_eq!(Subtype::try_from(b'S'), Ok(Subtype::UInt16));
        assert_eq!(Subtype::try_from(b'i'), Ok(Subtype::Int32));
        assert_eq!(Subtype::try_from(b'I'), Ok(Subtype::UInt32));
        assert_eq!(Subtype::try_from(b'f'), Ok(Subtype::Float));
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
