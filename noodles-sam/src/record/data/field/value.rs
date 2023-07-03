//! SAM record data field value and types.

pub mod array;
mod base_modifications;
pub mod character;
pub mod hex;

pub use self::{array::Array, character::Character, hex::Hex};

use std::{
    error,
    fmt::{self, Write},
    io, num,
};

use super::Type;

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
    Hex(Hex),
    /// An array (`B`).
    Array(Array),
}

impl Value {
    /// Parses a raw value as the given type.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam::record::data::field::{Type, Value};
    /// let value = Value::from_str_type("rg0", Type::String)?;
    /// assert_eq!(value, Value::String(String::from("rg0")));
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn from_str_type(s: &str, ty: Type) -> io::Result<Self> {
        use crate::reader::record::data::field::parse_value;

        let mut src = s.as_bytes();
        parse_value(&mut src, ty).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Returns the type of the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{Type, Value};
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
            Self::Array(_) => Type::Array,
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
    /// use noodles_sam::record::data::field::{value::Hex, Value};
    /// let hex: Hex = "CAFE".parse()?;
    /// assert_eq!(Value::Hex(hex.clone()).as_hex(), Some(&hex));
    /// assert_eq!(Value::Int32(0).as_hex(), None);
    /// # Ok::<_, noodles_sam::record::data::field::value::hex::ParseError>(())
    /// ```
    pub fn as_hex(&self) -> Option<&Hex> {
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
    /// assert!(Value::Hex("CAFE".parse()?).is_hex());
    /// assert!(!Value::Int32(0).is_hex());
    /// # Ok::<_, noodles_sam::record::data::field::value::hex::ParseError>(())
    /// ```
    pub fn is_hex(&self) -> bool {
        matches!(self, Self::Hex(_))
    }

    /// Returns the value as an array of 8-bit integers if it is an array of 8-bit integers.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Array, Value};
    /// let array = Array::UInt8(vec![0]);
    /// assert_eq!(Value::Array(array.clone()).as_array(), Some(&array));
    /// assert_eq!(Value::Int32(0).as_array(), None);
    /// ```
    pub fn as_array(&self) -> Option<&Array> {
        match *self {
            Self::Array(ref array) => Some(array),
            _ => None,
        }
    }

    /// Returns whether the value is an 8-bit integer array.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Array, Value};
    /// assert!(Value::Array(Array::UInt8(vec![0])).is_array());
    /// assert!(!Value::Int32(0).is_array());
    /// ```
    pub fn is_array(&self) -> bool {
        matches!(self, Self::Array(_))
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
            Self::Int8(n) => write!(f, "{n}"),
            Self::UInt8(n) => write!(f, "{n}"),
            Self::Int16(n) => write!(f, "{n}"),
            Self::UInt16(n) => write!(f, "{n}"),
            Self::Int32(n) => write!(f, "{n}"),
            Self::UInt32(n) => write!(f, "{n}"),
            Self::Float(n) => write!(f, "{n}"),
            Self::String(s) => f.write_str(s),
            Self::Hex(s) => write!(f, "{s}"),
            Self::Array(array) => write!(f, "{array}"),
        }
    }
}

/// An error returned when a raw SAM record data field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
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
    InvalidSubtype(array::subtype::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidIntValue(e) => Some(e),
            Self::InvalidFloatValue(e) => Some(e),
            Self::InvalidSubtype(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidCharacterValue => f.write_str("invalid character value"),
            Self::InvalidIntValue(_) => f.write_str("invalid int value"),
            Self::InvalidFloatValue(_) => f.write_str("invalid float value"),
            Self::InvalidStringValue => f.write_str("invalid string value"),
            Self::InvalidHexValue => f.write_str("invalid hex value"),
            Self::MissingSubtype => f.write_str("missing subtype"),
            Self::InvalidSubtype(_) => f.write_str("invalid subtype"),
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
    fn test_ty() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(
            Value::Character(Character::try_from('n')?).ty(),
            Type::Character
        );
        assert_eq!(Value::Int32(0).ty(), Type::Int32);
        assert_eq!(Value::Float(0.0).ty(), Type::Float);
        assert_eq!(Value::String(String::from("noodles")).ty(), Type::String);
        assert_eq!(Value::Hex("CAFE".parse()?).ty(), Type::Hex);
        assert_eq!(Value::Array(Array::UInt8(vec![0])).ty(), Type::Array);

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
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(Value::Character(Character::try_from('n')?).to_string(), "n");
        assert_eq!(Value::Int32(13).to_string(), "13");
        assert_eq!(Value::Float(0.0).to_string(), "0");

        assert_eq!(Value::String(String::new()).to_string(), "");
        assert_eq!(
            Value::String(String::from("noodles")).to_string(),
            "noodles"
        );

        assert_eq!(Value::Hex(Hex::default()).to_string(), "");
        assert_eq!(Value::Hex("CAFE".parse()?).to_string(), "CAFE");

        assert_eq!(
            Value::Array(Array::UInt8(vec![8, 13])).to_string(),
            "C,8,13"
        );

        Ok(())
    }
}
