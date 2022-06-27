//! SAM record data field value and types.

pub mod character;
pub mod subtype;
pub mod ty;

pub use self::{character::Character, subtype::Subtype, ty::Type};

use std::{
    error,
    fmt::{self, Write},
    io, num,
};

/// A SAM record data field value.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A character (`A`).
    Character(Character),
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
    String(String),
    /// A hex string (`H`).
    Hex(String),
    /// An 8-bit integer array (`Bc`).
    Int8Array(Vec<i8>),
    /// An 8-bit unsigned integer array (`BC`).
    UInt8Array(Vec<u8>),
    /// A 16-bit integer array (`Bs`).
    Int16Array(Vec<i16>),
    /// A 16-bit unsigned integer array (`BS`).
    UInt16Array(Vec<u16>),
    /// A 32-bit integer array (`Bi`).
    Int32Array(Vec<i32>),
    /// A 32-bit unsigned integer array (`BI`).
    UInt32Array(Vec<u32>),
    /// A single-precision floating-point array (`Bf`).
    FloatArray(Vec<f32>),
}

impl Value {
    /// Parses a raw value as the given type.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam::record::data::field::{value::Type, Value};
    /// let value = Value::from_str_type("rg0", Type::String)?;
    /// assert_eq!(value, Value::String(String::from("rg0")));
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn from_str_type(s: &str, ty: Type) -> io::Result<Self> {
        use crate::reader::record::data::field::parse_value;

        let mut src = s.as_bytes();
        parse_value(&mut src, ty)
    }

    /// Returns the type of the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Type, Value};
    /// assert_eq!(Value::Int32(0).ty(), Type::Int32);
    /// ```
    pub fn ty(&self) -> Type {
        match *self {
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
            Self::Int8Array(_) => Type::Array,
            Self::UInt8Array(_) => Type::Array,
            Self::Int16Array(_) => Type::Array,
            Self::UInt16Array(_) => Type::Array,
            Self::Int32Array(_) => Type::Array,
            Self::UInt32Array(_) => Type::Array,
            Self::FloatArray(_) => Type::Array,
        }
    }

    /// Returns the subtype of the value.
    ///
    /// Only arrays have subtypes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Subtype, Value};
    /// assert_eq!(Value::UInt8Array(vec![0]).subtype(), Some(Subtype::UInt8));
    /// assert_eq!(Value::Int32(0).subtype(), None);
    /// ```
    pub fn subtype(&self) -> Option<Subtype> {
        match *self {
            Self::Int8Array(_) => Some(Subtype::Int8),
            Self::UInt8Array(_) => Some(Subtype::UInt8),
            Self::Int16Array(_) => Some(Subtype::Int16),
            Self::UInt16Array(_) => Some(Subtype::UInt16),
            Self::Int32Array(_) => Some(Subtype::Int32),
            Self::UInt32Array(_) => Some(Subtype::UInt32),
            Self::FloatArray(_) => Some(Subtype::Float),
            _ => None,
        }
    }

    /// Returns the value as a character if it is a character.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Character, Value};
    /// let c = Character::try_from('a')?;
    /// assert_eq!(Value::Character(c).as_char(), Some('a'));
    /// assert_eq!(Value::Int32(0).as_char(), None);
    /// # Ok::<_, noodles_sam::record::data::field::value::character::ParseError>(())
    /// ```
    #[deprecated(since = "0.16.0", note = "Use `Value::as_character` instead.")]
    pub fn as_char(&self) -> Option<char> {
        match *self {
            Self::Character(c) => Some(char::from(c)),
            _ => None,
        }
    }

    /// Returns the value as a character if it is a character.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Character, Value};
    /// let c = Character::try_from('a')?;
    /// assert_eq!(Value::Character(c).as_character(), Some(c));
    /// assert_eq!(Value::Int32(0).as_character(), None);
    /// # Ok::<_, noodles_sam::record::data::field::value::character::ParseError>(())
    /// ```
    pub fn as_character(&self) -> Option<Character> {
        match *self {
            Self::Character(c) => Some(c),
            _ => None,
        }
    }

    /// Returns whether the value is a character.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Character, Value};
    /// assert!(Value::Character(Character::try_from('a')?).is_char());
    /// assert!(!Value::Int32(0).is_char());
    /// # Ok::<_, noodles_sam::record::data::field::value::character::ParseError>(())
    /// ```
    #[deprecated(since = "0.16.0", note = "Use `Value::is_character` instead.")]
    pub fn is_char(&self) -> bool {
        self.is_character()
    }

    /// Returns whether the value is a character.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Character, Value};
    /// assert!(Value::Character(Character::try_from('a')?).is_character());
    /// assert!(!Value::Int32(0).is_character());
    /// # Ok::<_, noodles_sam::record::data::field::value::character::ParseError>(())
    /// ```
    pub fn is_character(&self) -> bool {
        matches!(self, Self::Character(_))
    }

    /// Returns the value as a 32-bit integer if it is an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Int32(0).as_int32(), Some(0));
    /// assert_eq!(Value::Float(0.0).as_int32(), None);
    /// ```
    pub fn as_int32(&self) -> Option<i32> {
        match *self {
            Self::Int32(i) => Some(i),
            _ => None,
        }
    }

    /// Returns whether the value is an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Int32(0).is_int32());
    /// assert!(!Value::Float(0.0).is_int32());
    /// ```
    pub fn is_int32(&self) -> bool {
        matches!(self, Self::Int32(_))
    }

    /// Returns the value as a 64-bit integer.
    ///
    /// This is a convenience method that converts any integer to an `i64`, which captures the
    /// entire range of all record data field integer values.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Int32(0).as_int(), Some(0));
    /// assert_eq!(Value::Float(0.0).as_int(), None);
    /// ```
    pub fn as_int(&self) -> Option<i64> {
        match *self {
            Self::Int8(n) => Some(i64::from(n)),
            Self::UInt8(n) => Some(i64::from(n)),
            Self::Int16(n) => Some(i64::from(n)),
            Self::UInt16(n) => Some(i64::from(n)),
            Self::Int32(n) => Some(i64::from(n)),
            Self::UInt32(n) => Some(i64::from(n)),
            _ => None,
        }
    }

    /// Returns whether the value is an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
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

    /// Returns the value as a single-precision floating-point if it is a single-precision
    /// float-point.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Float(0.0).as_float(), Some(0.0));
    /// assert_eq!(Value::Int32(0).as_float(), None);
    /// ```
    pub fn as_float(&self) -> Option<f32> {
        match *self {
            Self::Float(s) => Some(s),
            _ => None,
        }
    }

    /// Returns whether the value is a single-precision floating-point.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Float(0.0).is_float());
    /// assert!(!Value::Int32(0).is_float());
    /// ```
    pub fn is_float(&self) -> bool {
        matches!(self, Self::Float(_))
    }

    /// Returns the value as a string slice if it is a string slice.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::String(String::from("noodles")).as_str(), Some("noodles"));
    /// assert_eq!(Value::Int32(0).as_str(), None);
    /// ```
    pub fn as_str(&self) -> Option<&str> {
        match *self {
            Self::String(ref s) => Some(s),
            _ => None,
        }
    }

    /// Returns whether the value is a string.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::String(String::from("noodles")).is_str());
    /// assert!(!Value::Int32(0).is_str());
    /// ```
    pub fn is_str(&self) -> bool {
        matches!(self, Self::String(_))
    }

    /// Returns the value as a string slice of hex if it is a string slice of hex.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Hex(String::from("CAFE")).as_hex(), Some("CAFE"));
    /// assert_eq!(Value::Int32(0).as_hex(), None);
    /// ```
    pub fn as_hex(&self) -> Option<&str> {
        match *self {
            Self::Hex(ref h) => Some(h),
            _ => None,
        }
    }

    /// Returns whether the value is a hex string.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Hex(String::from("CAFE")).is_hex());
    /// assert!(!Value::Int32(0).is_hex());
    /// ```
    pub fn is_hex(&self) -> bool {
        matches!(self, Self::Hex(_))
    }

    /// Returns the value as an array of 8-bit integers if it is an array of 8-bit integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Int8Array(vec![0]).as_int8_array(), Some(&[0][..]));
    /// assert_eq!(Value::Int32(0).as_int8_array(), None);
    /// ```
    pub fn as_int8_array(&self) -> Option<&[i8]> {
        match *self {
            Self::Int8Array(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is an 8-bit integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Int8Array(vec![0]).is_int8_array());
    /// assert!(!Value::Int32(0).is_int8_array());
    /// ```
    pub fn is_int8_array(&self) -> bool {
        matches!(self, Self::Int8Array(_))
    }

    /// Returns the value as an array of 8-bit unsigned integers if it is an array of 8-bit
    /// unsigned integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::UInt8Array(vec![0]).as_uint8_array(), Some(&[0][..]));
    /// assert_eq!(Value::Int32(0).as_uint8_array(), None);
    /// ```
    pub fn as_uint8_array(&self) -> Option<&[u8]> {
        match *self {
            Self::UInt8Array(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is an 8-bit unsigned integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::UInt8Array(vec![0]).is_uint8_array());
    /// assert!(!Value::Int32(0).is_uint8_array());
    /// ```
    pub fn is_uint8_array(&self) -> bool {
        matches!(self, Self::UInt8Array(_))
    }

    /// Returns the value as an array of 16-bit integers if it is an array of 16-bit integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Int16Array(vec![0]).as_int16_array(), Some(&[0][..]));
    /// assert_eq!(Value::Int32(0).as_int16_array(), None);
    /// ```
    pub fn as_int16_array(&self) -> Option<&[i16]> {
        match *self {
            Self::Int16Array(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is a 16-bit integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Int16Array(vec![0]).is_int16_array());
    /// assert!(!Value::Int32(0).is_int16_array());
    /// ```
    pub fn is_int16_array(&self) -> bool {
        matches!(self, Self::Int16Array(_))
    }

    /// Returns the value as an array of 16-bit unsigned integers if it is an array of 16-bit
    /// unsigned integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::UInt16Array(vec![0]).as_uint16_array(), Some(&[0][..]));
    /// assert_eq!(Value::Int32(0).as_uint16_array(), None);
    /// ```
    pub fn as_uint16_array(&self) -> Option<&[u16]> {
        match *self {
            Self::UInt16Array(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is a 16-bit unsigned integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::UInt16Array(vec![0]).is_uint16_array());
    /// assert!(!Value::Int32(0).is_int16_array());
    /// ```
    pub fn is_uint16_array(&self) -> bool {
        matches!(self, Self::UInt16Array(_))
    }

    /// Returns the value as an array of 32-bit integers if it is an array of 32-bit integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Int32Array(vec![0]).as_int32_array(), Some(&[0][..]));
    /// assert_eq!(Value::Int32(0).as_int32_array(), None);
    /// ```
    pub fn as_int32_array(&self) -> Option<&[i32]> {
        match *self {
            Self::Int32Array(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is a 32-bit integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Int32Array(vec![0]).is_int32_array());
    /// assert!(!Value::Int32(0).is_int16_array());
    /// ```
    pub fn is_int32_array(&self) -> bool {
        matches!(self, Self::Int32Array(_))
    }

    /// Returns the value as an array of 32-bit unsigned integers if it is an array of 32-bit
    /// unsigned integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::UInt32Array(vec![0]).as_uint32_array(), Some(&[0][..]));
    /// assert_eq!(Value::Int32(0).as_uint32_array(), None);
    /// ```
    pub fn as_uint32_array(&self) -> Option<&[u32]> {
        match *self {
            Self::UInt32Array(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is a 32-bit unsigned integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::UInt32Array(vec![0]).is_uint32_array());
    /// assert!(!Value::Int32(0).is_int32_array());
    /// ```
    pub fn is_uint32_array(&self) -> bool {
        matches!(self, Self::UInt32Array(_))
    }

    /// Returns the value as an array of single-precision floating-points if it is an array of
    /// single-precision floating-points.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::FloatArray(vec![0.0]).as_float_array(), Some(&[0.0][..]));
    /// assert_eq!(Value::Int32(0).as_float_array(), None);
    /// ```
    pub fn as_float_array(&self) -> Option<&[f32]> {
        match *self {
            Self::FloatArray(ref a) => Some(a),
            _ => None,
        }
    }

    /// Returns whether the value is a 32-bit integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Int32Array(vec![0]).is_int32_array());
    /// assert!(!Value::Int32(0).is_int16_array());
    /// ```
    pub fn is_float_array(&self) -> bool {
        matches!(self, Self::FloatArray(_))
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

impl From<Vec<i8>> for Value {
    fn from(values: Vec<i8>) -> Self {
        Value::Int8Array(values)
    }
}

impl From<Vec<u8>> for Value {
    fn from(values: Vec<u8>) -> Self {
        Value::UInt8Array(values)
    }
}

impl From<Vec<i16>> for Value {
    fn from(values: Vec<i16>) -> Self {
        Value::Int16Array(values)
    }
}

impl From<Vec<u16>> for Value {
    fn from(values: Vec<u16>) -> Self {
        Value::UInt16Array(values)
    }
}

impl From<Vec<i32>> for Value {
    fn from(values: Vec<i32>) -> Self {
        Value::Int32Array(values)
    }
}

impl From<Vec<u32>> for Value {
    fn from(values: Vec<u32>) -> Self {
        Value::UInt32Array(values)
    }
}

impl From<Vec<f32>> for Value {
    fn from(values: Vec<f32>) -> Self {
        Value::FloatArray(values)
    }
}

impl TryFrom<char> for Value {
    type Error = ParseError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        Character::try_from(c)
            .map(Self::Character)
            .map_err(|_| ParseError::InvalidCharacterValue)
    }
}

impl TryFrom<String> for Value {
    type Error = ParseError;

    fn try_from(s: String) -> Result<Self, Self::Error> {
        if is_valid_string(&s) {
            Ok(Self::String(s))
        } else {
            Err(ParseError::InvalidStringValue)
        }
    }
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Character(c) => f.write_char(char::from(*c)),
            Self::Int8(n) => write!(f, "{}", n),
            Self::UInt8(n) => write!(f, "{}", n),
            Self::Int16(n) => write!(f, "{}", n),
            Self::UInt16(n) => write!(f, "{}", n),
            Self::Int32(n) => write!(f, "{}", n),
            Self::UInt32(n) => write!(f, "{}", n),
            Self::Float(n) => write!(f, "{}", n),
            Self::String(s) => f.write_str(s),
            Self::Hex(s) => f.write_str(s),
            Self::Int8Array(values) => {
                f.write_char(char::from(Subtype::Int8))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
            Self::UInt8Array(values) => {
                f.write_char(char::from(Subtype::UInt8))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
            Self::Int16Array(values) => {
                f.write_char(char::from(Subtype::Int16))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
            Self::UInt16Array(values) => {
                f.write_char(char::from(Subtype::UInt16))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
            Self::Int32Array(values) => {
                f.write_char(char::from(Subtype::Int32))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
            Self::UInt32Array(values) => {
                f.write_char(char::from(Subtype::UInt32))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
            Self::FloatArray(values) => {
                f.write_char(char::from(Subtype::Float))?;

                for value in values {
                    write!(f, ",{}", value)?;
                }

                Ok(())
            }
        }
    }
}

/// An error returned when a raw SAM record data field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The data field type is unsupported.
    UnsupportedType(Type),
    /// The data field character value is invalid.
    InvalidCharacterValue,
    /// The data field integer value is invalid.
    InvalidIntValue(num::ParseIntError),
    /// The data field floating-point value is invalid.
    InvalidFloatValue(num::ParseFloatError),
    /// The data field string value is invalid.
    InvalidStringValue,
    /// The data field hex value is invalid.
    InvalidHexValue,
    /// The data field subtype is missing.
    MissingSubtype,
    /// The data field subtype is invalid.
    InvalidSubtype(subtype::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::UnsupportedType(ty) => write!(f, "unsupported type: {}", ty),
            Self::InvalidCharacterValue => f.write_str("invalid character value"),
            Self::InvalidIntValue(e) => write!(f, "invalid int value: {}", e),
            Self::InvalidFloatValue(e) => write!(f, "invalid float value: {}", e),
            Self::InvalidStringValue => write!(f, "invalid string value"),
            Self::InvalidHexValue => write!(f, "invalid hex value"),
            Self::MissingSubtype => f.write_str("missing subtype"),
            Self::InvalidSubtype(e) => write!(f, "invalid subtype: {}", e),
        }
    }
}

// § 1.5 The alignment section: optional fields (2021-01-07)
fn is_valid_string_char(c: char) -> bool {
    matches!(c, ' ' | '!'..='~')
}

fn is_valid_string(s: &str) -> bool {
    s.chars().all(is_valid_string_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ty() -> Result<(), character::ParseError> {
        assert_eq!(
            Value::Character(Character::try_from('n')?).ty(),
            Type::Character
        );
        assert_eq!(Value::Int32(0).ty(), Type::Int32);
        assert_eq!(Value::Float(0.0).ty(), Type::Float);
        assert_eq!(Value::String(String::from("noodles")).ty(), Type::String);
        assert_eq!(Value::Hex(String::from("CAFE")).ty(), Type::Hex);
        assert_eq!(Value::Int8Array(vec![0]).ty(), Type::Array);
        assert_eq!(Value::UInt8Array(vec![0]).ty(), Type::Array);
        assert_eq!(Value::Int16Array(vec![0]).ty(), Type::Array);
        assert_eq!(Value::UInt16Array(vec![0]).ty(), Type::Array);
        assert_eq!(Value::Int32Array(vec![0]).ty(), Type::Array);
        assert_eq!(Value::UInt32Array(vec![0]).ty(), Type::Array);
        assert_eq!(Value::FloatArray(vec![0.0]).ty(), Type::Array);

        Ok(())
    }

    #[test]
    fn test_subtype() -> Result<(), character::ParseError> {
        assert_eq!(Value::Character(Character::try_from('n')?).subtype(), None);
        assert_eq!(Value::Int32(0).subtype(), None);
        assert_eq!(Value::Float(0.0).subtype(), None);
        assert_eq!(Value::String(String::from("noodles")).subtype(), None);
        assert_eq!(Value::Hex(String::from("CAFE")).subtype(), None);
        assert_eq!(Value::Int8Array(vec![0]).subtype(), Some(Subtype::Int8));
        assert_eq!(Value::UInt8Array(vec![0]).subtype(), Some(Subtype::UInt8));
        assert_eq!(Value::Int16Array(vec![0]).subtype(), Some(Subtype::Int16));
        assert_eq!(Value::UInt16Array(vec![0]).subtype(), Some(Subtype::UInt16));
        assert_eq!(Value::Int32Array(vec![0]).subtype(), Some(Subtype::Int32));
        assert_eq!(Value::UInt32Array(vec![0]).subtype(), Some(Subtype::UInt32));
        assert_eq!(Value::FloatArray(vec![0.0]).subtype(), Some(Subtype::Float));

        Ok(())
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
        assert_eq!(Value::from(vec![0i8]), Value::Int8Array(vec![0]));
    }

    #[test]
    fn test_from_vec_u8_for_value() {
        assert_eq!(Value::from(vec![0u8]), Value::UInt8Array(vec![0]));
    }

    #[test]
    fn test_from_vec_i16_for_value() {
        assert_eq!(Value::from(vec![0i16]), Value::Int16Array(vec![0]));
    }

    #[test]
    fn test_from_vec_u16_for_value() {
        assert_eq!(Value::from(vec![0u16]), Value::UInt16Array(vec![0]));
    }

    #[test]
    fn test_from_vec_i32_for_value() {
        assert_eq!(Value::from(vec![0i32]), Value::Int32Array(vec![0]));
    }

    #[test]
    fn test_from_vec_u32_for_value() {
        assert_eq!(Value::from(vec![0u32]), Value::UInt32Array(vec![0]));
    }

    #[test]
    fn test_from_vec_f32_for_value() {
        assert_eq!(Value::from(vec![0.0f32]), Value::FloatArray(vec![0.0]));
    }

    #[test]
    fn test_try_from_char_for_value() -> Result<(), character::ParseError> {
        assert_eq!(
            Value::try_from('n'),
            Ok(Value::Character(Character::try_from('n')?))
        );

        assert_eq!(
            Value::try_from('🍜'),
            Err(ParseError::InvalidCharacterValue)
        );

        Ok(())
    }

    #[test]
    fn test_try_from_string_for_value() {
        assert_eq!(
            Value::try_from(String::from("noodles")),
            Ok(Value::String(String::from("noodles")))
        );
    }

    #[test]
    fn test_fmt() -> Result<(), character::ParseError> {
        assert_eq!(Value::Character(Character::try_from('n')?).to_string(), "n");
        assert_eq!(Value::Int32(13).to_string(), "13");
        assert_eq!(Value::Float(0.0).to_string(), "0");

        assert_eq!(Value::String(String::new()).to_string(), "");
        assert_eq!(
            Value::String(String::from("noodles")).to_string(),
            "noodles"
        );

        assert_eq!(Value::Hex(String::new()).to_string(), "");
        assert_eq!(Value::Hex(String::from("CAFE")).to_string(), "CAFE");

        assert_eq!(Value::Int8Array(Vec::new()).to_string(), "c");
        assert_eq!(Value::Int8Array(vec![1, -2]).to_string(), "c,1,-2");

        assert_eq!(Value::UInt8Array(Vec::new()).to_string(), "C");
        assert_eq!(Value::UInt8Array(vec![3, 5]).to_string(), "C,3,5");

        assert_eq!(Value::Int16Array(Vec::new()).to_string(), "s");
        assert_eq!(Value::Int16Array(vec![8, -13]).to_string(), "s,8,-13");

        assert_eq!(Value::UInt16Array(Vec::new()).to_string(), "S");
        assert_eq!(Value::UInt16Array(vec![21, 34]).to_string(), "S,21,34");

        assert_eq!(Value::Int32Array(Vec::new()).to_string(), "i");
        assert_eq!(Value::Int32Array(vec![55, -89]).to_string(), "i,55,-89");

        assert_eq!(Value::UInt32Array(Vec::new()).to_string(), "I");
        assert_eq!(Value::UInt32Array(vec![144, 233]).to_string(), "I,144,233");

        assert_eq!(Value::FloatArray(Vec::new()).to_string(), "f");
        assert_eq!(Value::FloatArray(vec![0.0, 1.0]).to_string(), "f,0,1");

        Ok(())
    }
}
