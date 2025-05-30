//! Alignment record data field value buffer.

mod array;

use bstr::BString;

pub use self::array::Array;
use crate::alignment::record::data::field::Type;

/// An alignment record data field value buffer.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A character (`A`).
    Character(u8),
    /// An 8-bit integer (`c`).
    Int8(i8),
    /// An 8-bit unsigned integer (`C`).
    UInt8(u8),
    /// A 16-bit integer (`s`).
    Int16(i16),
    /// A 16-bit unsigned integer (`S`).
    UInt16(u16),
    /// A 32-bit integer (`i`).
    Int32(i32),
    /// A 32-bit unsigned integer (`I`).
    UInt32(u32),
    /// A single-precision floating-point (`f`).
    Float(f32),
    /// A string (`Z`).
    String(BString),
    /// A hex string (`H`).
    Hex(BString),
    /// An array (`B`).
    Array(Array),
}

impl Value {
    /// Returns the type of the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::{record::data::field::Type, record_buf::data::field::Value};
    /// assert_eq!(Value::Int32(0).ty(), Type::Int32);
    /// ```
    pub fn ty(&self) -> Type {
        match self {
            Self::Character(_) => Type::Character,
            Self::Int8(_) => Type::Int8,
            Self::UInt8(_) => Type::UInt8,
            Self::Int16(_) => Type::Int16,
            Self::UInt16(_) => Type::UInt16,
            Self::Int32(_) => Type::Int32,
            Self::UInt32(_) => Type::UInt32,
            Self::Float(_) => Type::Float,
            Self::String(_) => Type::String,
            Self::Hex(_) => Type::Hex,
            Self::Array(_) => Type::Array,
        }
    }

    /// Returns the value as a 64-bit integer.
    ///
    /// This is a convenience method that converts any integer to an `i64`, which captures the
    /// entire range of all alignment record data field integer values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record_buf::data::field::Value;
    /// assert_eq!(Value::from(0).as_int(), Some(0));
    /// assert_eq!(Value::from("noodles").as_int(), None);
    /// ```
    pub fn as_int(&self) -> Option<i64> {
        match self {
            Self::Int8(n) => Some(i64::from(*n)),
            Self::UInt8(n) => Some(i64::from(*n)),
            Self::Int16(n) => Some(i64::from(*n)),
            Self::UInt16(n) => Some(i64::from(*n)),
            Self::Int32(n) => Some(i64::from(*n)),
            Self::UInt32(n) => Some(i64::from(*n)),
            _ => None,
        }
    }

    /// Returns whether the value is an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record_buf::data::field::Value;
    /// assert!(Value::Int32(0).is_int());
    /// assert!(!Value::Float(0.0).is_int());
    /// ```
    pub fn is_int(&self) -> bool {
        matches!(
            self,
            Self::Int8(_)
                | Self::UInt8(_)
                | Self::Int16(_)
                | Self::UInt16(_)
                | Self::Int32(_)
                | Self::UInt32(_)
        )
    }
}

impl From<i8> for Value {
    fn from(n: i8) -> Self {
        if n >= 0 {
            Self::from(n as u8)
        } else {
            Self::Int8(n)
        }
    }
}

impl From<u8> for Value {
    fn from(n: u8) -> Self {
        Self::UInt8(n)
    }
}

impl From<i16> for Value {
    fn from(n: i16) -> Self {
        if n >= 0 {
            Self::from(n as u16)
        } else if n >= i16::from(i8::MIN) {
            Self::Int8(n as i8)
        } else {
            Self::Int16(n)
        }
    }
}

impl From<u16> for Value {
    fn from(n: u16) -> Self {
        if n <= u16::from(u8::MAX) {
            Self::UInt8(n as u8)
        } else {
            Self::UInt16(n)
        }
    }
}

impl From<i32> for Value {
    fn from(n: i32) -> Self {
        if n >= 0 {
            Self::from(n as u32)
        } else if n >= i32::from(i8::MIN) {
            Self::Int8(n as i8)
        } else if n >= i32::from(i16::MIN) {
            Self::Int16(n as i16)
        } else {
            Self::Int32(n)
        }
    }
}

impl From<u32> for Value {
    fn from(n: u32) -> Self {
        if n <= u32::from(u8::MAX) {
            Self::UInt8(n as u8)
        } else if n <= u32::from(u16::MAX) {
            Self::UInt16(n as u16)
        } else {
            Self::UInt32(n)
        }
    }
}

impl From<f32> for Value {
    fn from(n: f32) -> Self {
        Value::Float(n)
    }
}

impl From<&str> for Value {
    fn from(s: &str) -> Self {
        Self::String(s.into())
    }
}

impl From<String> for Value {
    fn from(s: String) -> Self {
        Self::String(s.into())
    }
}

impl From<Vec<i8>> for Value {
    fn from(values: Vec<i8>) -> Self {
        Value::Array(Array::Int8(values))
    }
}

impl From<Vec<u8>> for Value {
    fn from(values: Vec<u8>) -> Self {
        Value::Array(Array::UInt8(values))
    }
}

impl From<Vec<i16>> for Value {
    fn from(values: Vec<i16>) -> Self {
        Value::Array(Array::Int16(values))
    }
}

impl From<Vec<u16>> for Value {
    fn from(values: Vec<u16>) -> Self {
        Value::Array(Array::UInt16(values))
    }
}

impl From<Vec<i32>> for Value {
    fn from(values: Vec<i32>) -> Self {
        Value::Array(Array::Int32(values))
    }
}

impl From<Vec<u32>> for Value {
    fn from(values: Vec<u32>) -> Self {
        Value::Array(Array::UInt32(values))
    }
}

impl From<Vec<f32>> for Value {
    fn from(values: Vec<f32>) -> Self {
        Value::Array(Array::Float(values))
    }
}

impl TryFrom<i64> for Value {
    type Error = crate::io::reader::record_buf::data::field::value::ParseError;

    fn try_from(n: i64) -> Result<Self, Self::Error> {
        const MIN: i64 = i32::MIN as i64;
        const MAX: i64 = u32::MAX as i64;

        if n > MAX {
            Err(Self::Error::InvalidIntegerValue)
        } else if n >= 0 {
            Ok(Self::from(n as u32))
        } else if n >= i64::from(i8::MIN) {
            Ok(Self::Int8(n as i8))
        } else if n >= i64::from(i16::MIN) {
            Ok(Self::Int16(n as i16))
        } else if n >= MIN {
            Ok(Self::Int32(n as i32))
        } else {
            Err(Self::Error::InvalidIntegerValue)
        }
    }
}

impl<'a> From<&'a Value> for crate::alignment::record::data::field::Value<'a> {
    fn from(value_buf: &'a Value) -> Self {
        match value_buf {
            Value::Character(c) => Self::Character(*c),
            Value::Int8(n) => Self::Int8(*n),
            Value::UInt8(n) => Self::UInt8(*n),
            Value::Int16(n) => Self::Int16(*n),
            Value::UInt16(n) => Self::UInt16(*n),
            Value::Int32(n) => Self::Int32(*n),
            Value::UInt32(n) => Self::UInt32(*n),
            Value::Float(n) => Self::Float(*n),
            Value::String(s) => Self::String(s.as_ref()),
            Value::Hex(s) => Self::Hex(s.as_ref()),
            Value::Array(array) => Self::Array(array.into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ty() {
        assert_eq!(Value::Character(b'n').ty(), Type::Character);
        assert_eq!(Value::Int32(0).ty(), Type::Int32);
        assert_eq!(Value::Float(0.0).ty(), Type::Float);
        assert_eq!(Value::from("noodles").ty(), Type::String);
        assert_eq!(Value::Hex(b"CAFE".into()).ty(), Type::Hex);
        assert_eq!(Value::Array(Array::UInt8(vec![0])).ty(), Type::Array);
    }

    #[test]
    fn test_from_i8_for_value() {
        assert_eq!(Value::from(i8::MIN), Value::Int8(i8::MIN));
        assert_eq!(Value::from(-1i8), Value::Int8(-1));

        assert_eq!(Value::from(0i8), Value::UInt8(0));
        assert_eq!(Value::from(1i8), Value::UInt8(1));
        assert_eq!(Value::from(i8::MAX), Value::UInt8(i8::MAX as u8));
    }

    #[test]
    fn test_from_u8_for_value() {
        assert_eq!(Value::from(u8::MIN), Value::UInt8(u8::MIN));
        assert_eq!(Value::from(u8::MAX), Value::UInt8(u8::MAX));
    }

    #[test]
    fn test_from_i16_for_value() {
        assert_eq!(Value::from(i16::MIN), Value::Int16(i16::MIN));
        assert_eq!(Value::from(-129i16), Value::Int16(-129)); // i8::MIN - 1

        assert_eq!(Value::from(-128i16), Value::Int8(-128)); // i8::MAX
        assert_eq!(Value::from(-1i16), Value::Int8(-1));

        assert_eq!(Value::from(0i16), Value::UInt8(0));
        assert_eq!(Value::from(1i16), Value::UInt8(1));
        assert_eq!(Value::from(255i16), Value::UInt8(255)); // u8::MAX

        assert_eq!(Value::from(256i16), Value::UInt16(256)); // u8::MAX + 1
        assert_eq!(Value::from(i16::MAX), Value::UInt16(i16::MAX as u16));
    }

    #[test]
    fn test_from_u16_for_value() {
        assert_eq!(Value::from(u16::MIN), Value::UInt8(0));
        assert_eq!(Value::from(255u16), Value::UInt8(255));

        assert_eq!(Value::from(256u16), Value::UInt16(256));
        assert_eq!(Value::from(u16::MAX), Value::UInt16(u16::MAX));
    }

    #[test]
    fn test_from_i32_for_value() {
        assert_eq!(Value::from(i32::MIN), Value::Int32(i32::MIN));
        assert_eq!(Value::from(-32769i32), Value::Int32(-32769)); // i16::MIN - 1

        assert_eq!(Value::from(-32768i32), Value::Int16(-32768)); // i16::MIN
        assert_eq!(Value::from(-129i32), Value::Int16(-129)); // i8::MIN - 1

        assert_eq!(Value::from(-128i32), Value::Int8(-128)); // i8::MIN
        assert_eq!(Value::from(-1i32), Value::Int8(-1));

        assert_eq!(Value::from(0i32), Value::UInt8(0));
        assert_eq!(Value::from(1i32), Value::UInt8(1));
        assert_eq!(Value::from(255i32), Value::UInt8(255)); // u8::MAX

        assert_eq!(Value::from(256i32), Value::UInt16(256)); // u8::MAX + 1
        assert_eq!(Value::from(65535i32), Value::UInt16(65535)); // u16::MAX

        assert_eq!(Value::from(65536i32), Value::UInt32(65536)); // u16::MAX + 1
        assert_eq!(Value::from(i32::MAX), Value::UInt32(i32::MAX as u32));
    }

    #[test]
    fn test_from_u32_for_value() {
        assert_eq!(Value::from(u32::MIN), Value::UInt8(0));
        assert_eq!(Value::from(255u32), Value::UInt8(255)); // u8::MAX

        assert_eq!(Value::from(256u32), Value::UInt16(256)); // u8::MAX + 1
        assert_eq!(Value::from(65535u32), Value::UInt16(65535)); // u16::MAX

        assert_eq!(Value::from(65536u32), Value::UInt32(65536)); // u16::MAX + 1
        assert_eq!(Value::from(u32::MAX), Value::UInt32(u32::MAX));
    }

    #[test]
    fn test_from_f32_for_value() {
        assert_eq!(Value::from(0.0f32), Value::Float(0.0));
    }

    #[test]
    fn test_from_vec_i8_for_value() {
        assert_eq!(Value::from(vec![0i8]), Value::Array(Array::Int8(vec![0])));
    }

    #[test]
    fn test_from_vec_u8_for_value() {
        assert_eq!(Value::from(vec![0u8]), Value::Array(Array::UInt8(vec![0])));
    }

    #[test]
    fn test_from_vec_i16_for_value() {
        assert_eq!(Value::from(vec![0i16]), Value::Array(Array::Int16(vec![0])));
    }

    #[test]
    fn test_from_vec_u16_for_value() {
        assert_eq!(
            Value::from(vec![0u16]),
            Value::Array(Array::UInt16(vec![0]))
        );
    }

    #[test]
    fn test_from_vec_i32_for_value() {
        assert_eq!(Value::from(vec![0i32]), Value::Array(Array::Int32(vec![0])));
    }

    #[test]
    fn test_from_vec_u32_for_value() {
        assert_eq!(
            Value::from(vec![0u32]),
            Value::Array(Array::UInt32(vec![0]))
        );
    }

    #[test]
    fn test_from_vec_f32_for_value() {
        assert_eq!(
            Value::from(vec![0.0f32]),
            Value::Array(Array::Float(vec![0.0]))
        );
    }

    #[test]
    fn test_try_from_i64_for_value()
    -> Result<(), crate::io::reader::record_buf::data::field::value::ParseError> {
        use crate::io::reader::record_buf::data::field::value::ParseError;

        fn t(n: i64, expected: Value) -> Result<(), ParseError> {
            let actual = Value::try_from(n)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        assert_eq!(
            Value::try_from(-2147483649i64),
            Err(ParseError::InvalidIntegerValue)
        );

        t(-2147483648, Value::Int32(i32::MIN))?;
        t(-2147483647, Value::Int32(-2147483647))?;

        t(-32769, Value::Int32(-32769))?;
        t(-32768, Value::Int16(i16::MIN))?;
        t(-32767, Value::Int16(-32767))?;

        t(-129, Value::Int16(-129))?;
        t(-128, Value::Int8(i8::MIN))?;
        t(-127, Value::Int8(-127))?;

        t(-1, Value::Int8(-1))?;
        t(0, Value::UInt8(0))?;
        t(1, Value::UInt8(1))?;

        t(254, Value::UInt8(254))?;
        t(255, Value::UInt8(u8::MAX))?;
        t(256, Value::UInt16(256))?;

        t(65534, Value::UInt16(65534))?;
        t(65535, Value::UInt16(u16::MAX))?;
        t(65536, Value::UInt32(65536))?;

        t(4294967294, Value::UInt32(4294967294))?;
        t(4294967295, Value::UInt32(u32::MAX))?;

        assert_eq!(
            Value::try_from(4294967296i64),
            Err(ParseError::InvalidIntegerValue)
        );

        Ok(())
    }
}
