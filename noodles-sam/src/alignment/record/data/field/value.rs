//! Alignment record data field value.

pub mod array;

pub use self::array::Array;
use super::Type;

/// An alignment record data field value.
#[derive(Debug)]
pub enum Value<'a> {
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
    String(&'a [u8]),
    /// A hex string (`H`).
    Hex(&'a [u8]),
    /// An array (`B`).
    Array(Array<'a>),
}

impl<'a> Value<'a> {
    /// Return the type of the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::data::field::{Type, Value};
    /// assert_eq!(Value::UInt8(0).ty(), Type::UInt8);
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
    /// entire range of all record data field integer values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::data::field::Value;
    /// assert_eq!(Value::Int8(8).as_int(), Some(8));
    /// assert_eq!(Value::UInt32(13).as_int(), Some(13));
    /// assert!(Value::Float(0.0).as_int().is_none());
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
}
