//! SAM record data field value and types.

pub mod subtype;
pub mod ty;

pub use self::{subtype::Subtype, ty::Type};

use std::{
    error,
    fmt::{self, Write},
    num,
    str::FromStr,
};

use super::DELIMITER;

const ARRAY_VALUE_DELIMITER: char = ',';

/// A SAM record data field value.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A character (`A`).
    Char(char),
    /// An integer (`i`).
    Int(i64),
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
    /// Returns the type of the value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::{value::Type, Value};
    /// assert_eq!(Value::Int(0).ty(), Type::Int);
    /// ```
    pub fn ty(&self) -> Type {
        match *self {
            Self::Char(_) => Type::Char,
            Self::Int(_) => Type::Int,
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
    /// assert_eq!(Value::Int(0).subtype(), None);
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
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Char('a').as_char(), Some('a'));
    /// assert_eq!(Value::Int(0).as_char(), None);
    /// ```
    pub fn as_char(&self) -> Option<char> {
        match *self {
            Self::Char(c) => Some(c),
            _ => None,
        }
    }

    /// Returns whether the value is a character.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Char('a').is_char());
    /// assert!(!Value::Int(0).is_char());
    /// ```
    pub fn is_char(&self) -> bool {
        matches!(self, Self::Char(_))
    }

    /// Returns the value as a 32-bit integer if it is an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Int(0).as_int(), Some(0));
    /// assert_eq!(Value::Char('a').as_int(), None);
    /// ```
    pub fn as_int(&self) -> Option<i64> {
        match *self {
            Self::Int(i) => Some(i),
            _ => None,
        }
    }

    /// Returns whether the value is an integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert!(Value::Int(0).is_int());
    /// assert!(!Value::Char('a').is_int());
    /// ```
    pub fn is_int(&self) -> bool {
        matches!(self, Self::Int(_))
    }

    /// Returns the value as a single-precision floating-point if it is a single-precision
    /// float-point.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::field::Value;
    /// assert_eq!(Value::Float(0.0).as_float(), Some(0.0));
    /// assert_eq!(Value::Int(0).as_float(), None);
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
    /// assert!(!Value::Int(0).is_float());
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
    /// assert_eq!(Value::Int(0).as_str(), None);
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
    /// assert!(!Value::Int(0).is_str());
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
    /// assert_eq!(Value::Int(0).as_hex(), None);
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
    /// assert!(!Value::Int(0).is_hex());
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
    /// assert_eq!(Value::Int(0).as_int8_array(), None);
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
    /// assert!(!Value::Int(0).is_int8_array());
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
    /// assert_eq!(Value::Int(0).as_uint8_array(), None);
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
    /// assert!(!Value::Int(0).is_uint8_array());
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
    /// assert_eq!(Value::Int(0).as_int16_array(), None);
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
    /// assert!(!Value::Int(0).is_int16_array());
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
    /// assert_eq!(Value::Int(0).as_uint16_array(), None);
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
    /// assert!(!Value::Int(0).is_int16_array());
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
    /// assert_eq!(Value::Int(0).as_int32_array(), None);
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
    /// assert!(!Value::Int(0).is_int16_array());
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
    /// assert_eq!(Value::Int(0).as_uint32_array(), None);
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
    /// assert!(!Value::Int(0).is_int32_array());
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
    /// assert_eq!(Value::Int(0).as_float_array(), None);
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
    /// assert!(!Value::Int(0).is_int16_array());
    /// ```
    pub fn is_float_array(&self) -> bool {
        matches!(self, Self::FloatArray(_))
    }
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Char(c) => f.write_char(*c),
            Self::Int(n) => write!(f, "{}", n),
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
    /// The data field type is invalid.
    InvalidType(ty::ParseError),
    /// The data field character value is invalid.
    InvalidCharValue,
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
            Self::InvalidType(e) => write!(f, "invalid type: {}", e),
            Self::InvalidCharValue => f.write_str("invalid char value"),
            Self::InvalidIntValue(e) => write!(f, "invalid int value: {}", e),
            Self::InvalidFloatValue(e) => write!(f, "invalid float value: {}", e),
            Self::InvalidStringValue => write!(f, "invalid string value"),
            Self::InvalidHexValue => write!(f, "invalid hex value"),
            Self::MissingSubtype => f.write_str("missing subtype"),
            Self::InvalidSubtype(e) => write!(f, "invalid subtype: {}", e),
        }
    }
}

impl FromStr for Value {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.split_once(DELIMITER) {
            Some((t, v)) => t
                .parse()
                .map_err(ParseError::InvalidType)
                .and_then(|ty| parse_value(ty, v)),
            None => Err(ParseError::Invalid),
        }
    }
}

fn parse_value(ty: Type, s: &str) -> Result<Value, ParseError> {
    match ty {
        Type::Char => parse_char(s).map(Value::Char),
        Type::Int => parse_int(s).map(Value::Int),
        Type::Float => parse_f32(s).map(Value::Float),
        Type::String => parse_string(s).map(Value::String),
        Type::Hex => parse_hex(s).map(Value::Hex),
        Type::Array => parse_array(s),
    }
}

fn parse_char(s: &str) -> Result<char, ParseError> {
    if let Some(c) = s.chars().next() {
        if c.is_ascii_graphic() {
            return Ok(c);
        }
    }

    Err(ParseError::InvalidCharValue)
}

fn parse_i8(s: &str) -> Result<i8, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_u8(s: &str) -> Result<u8, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_i16(s: &str) -> Result<i16, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_u16(s: &str) -> Result<u16, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_i32(s: &str) -> Result<i32, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_u32(s: &str) -> Result<u32, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_int(s: &str) -> Result<i64, ParseError> {
    s.parse().map_err(ParseError::InvalidIntValue)
}

fn parse_f32(s: &str) -> Result<f32, ParseError> {
    s.parse().map_err(ParseError::InvalidFloatValue)
}

// § 1.5 The alignment section: optional fields (2021-01-07)
fn is_valid_string_char(c: char) -> bool {
    matches!(c, ' ' | '!'..='~')
}

fn parse_string(s: &str) -> Result<String, ParseError> {
    if s.chars().all(is_valid_string_char) {
        Ok(s.into())
    } else {
        Err(ParseError::InvalidStringValue)
    }
}

// § 1.5 The alignment section: optional fields (2021-01-07)
fn is_valid_hex_char(c: char) -> bool {
    matches!(c, '0'..='9' | 'A'..='F')
}

fn parse_hex(s: &str) -> Result<String, ParseError> {
    if s.len() % 2 == 0 && s.chars().all(is_valid_hex_char) {
        Ok(s.into())
    } else {
        Err(ParseError::InvalidHexValue)
    }
}

fn parse_array(s: &str) -> Result<Value, ParseError> {
    let mut raw_values = s.split(ARRAY_VALUE_DELIMITER);

    let subtype = raw_values
        .next()
        .ok_or(ParseError::MissingSubtype)
        .and_then(|t| t.parse().map_err(ParseError::InvalidSubtype))?;

    match subtype {
        Subtype::Int8 => raw_values
            .map(parse_i8)
            .collect::<Result<_, _>>()
            .map(Value::Int8Array),

        Subtype::UInt8 => raw_values
            .map(parse_u8)
            .collect::<Result<_, _>>()
            .map(Value::UInt8Array),

        Subtype::Int16 => raw_values
            .map(parse_i16)
            .collect::<Result<_, _>>()
            .map(Value::Int16Array),

        Subtype::UInt16 => raw_values
            .map(parse_u16)
            .collect::<Result<_, _>>()
            .map(Value::UInt16Array),

        Subtype::Int32 => raw_values
            .map(parse_i32)
            .collect::<Result<_, _>>()
            .map(Value::Int32Array),

        Subtype::UInt32 => raw_values
            .map(parse_u32)
            .collect::<Result<_, _>>()
            .map(Value::UInt32Array),

        Subtype::Float => raw_values
            .map(parse_f32)
            .collect::<Result<_, _>>()
            .map(Value::FloatArray),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ty() {
        assert_eq!(Value::Char('n').ty(), Type::Char);
        assert_eq!(Value::Int(0).ty(), Type::Int);
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
    }

    #[test]
    fn test_subtype() {
        assert_eq!(Value::Char('n').subtype(), None);
        assert_eq!(Value::Int(0).subtype(), None);
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
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Value::Char('n').to_string(), "n");
        assert_eq!(Value::Int(13).to_string(), "13");
        assert_eq!(Value::Float(0.0).to_string(), "0");

        assert_eq!(
            Value::String(String::from("noodles")).to_string(),
            "noodles"
        );

        assert_eq!(Value::Hex(String::from("CAFE")).to_string(), "CAFE");
        assert_eq!(Value::Int8Array(vec![1, -2]).to_string(), "c,1,-2");
        assert_eq!(Value::UInt8Array(vec![3, 5]).to_string(), "C,3,5");
        assert_eq!(Value::Int16Array(vec![8, -13]).to_string(), "s,8,-13");
        assert_eq!(Value::UInt16Array(vec![21, 34]).to_string(), "S,21,34");
        assert_eq!(Value::Int32Array(vec![55, -89]).to_string(), "i,55,-89");
        assert_eq!(Value::UInt32Array(vec![144, 233]).to_string(), "I,144,233");
        assert_eq!(Value::FloatArray(vec![0.0, 1.0]).to_string(), "f,0,1");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("A:n".parse(), Ok(Value::Char('n')));
        assert_eq!("A:".parse::<Value>(), Err(ParseError::InvalidCharValue));
        assert_eq!("A:🍜".parse::<Value>(), Err(ParseError::InvalidCharValue));

        assert_eq!("i:13".parse(), Ok(Value::Int(13)));
        assert!(matches!(
            "i:".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "i:ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!("f:0.0".parse(), Ok(Value::Float(0.0)));
        assert!(matches!(
            "f:".parse::<Value>(),
            Err(ParseError::InvalidFloatValue(_))
        ));
        assert!(matches!(
            "f:ndls".parse::<Value>(),
            Err(ParseError::InvalidFloatValue(_))
        ));

        assert_eq!("Z:".parse(), Ok(Value::String(String::from(""))));
        assert_eq!("Z: ".parse(), Ok(Value::String(String::from(" "))));
        assert_eq!(
            "Z:noodles".parse(),
            Ok(Value::String(String::from("noodles")))
        );
        assert_eq!("Z:🍜".parse::<Value>(), Err(ParseError::InvalidStringValue));

        assert_eq!("H:CAFE".parse(), Ok(Value::Hex(String::from("CAFE"))));
        assert_eq!("H:cafe".parse::<Value>(), Err(ParseError::InvalidHexValue));
        assert_eq!("H:CAFE0".parse::<Value>(), Err(ParseError::InvalidHexValue));
        assert_eq!("H:NDLS".parse::<Value>(), Err(ParseError::InvalidHexValue));

        assert_eq!("B:c,1,-2".parse(), Ok(Value::Int8Array(vec![1, -2])));
        assert!(matches!(
            "B:c,".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "B:c,ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!("B:C,3,5".parse(), Ok(Value::UInt8Array(vec![3, 5])));
        assert!(matches!(
            "B:C,".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "B:C,ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!("B:s,8,-13".parse(), Ok(Value::Int16Array(vec![8, -13])));
        assert!(matches!(
            "B:s,".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "B:s,ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!("B:S,21,34".parse(), Ok(Value::UInt16Array(vec![21, 34])));
        assert!(matches!(
            "B:S,".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "B:S,ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!("B:i,55,-89".parse(), Ok(Value::Int32Array(vec![55, -89])));
        assert!(matches!(
            "B:i,".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "B:i,ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!(
            "B:I,144,233".parse(),
            Ok(Value::UInt32Array(vec![144, 233]))
        );
        assert!(matches!(
            "B:I,".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));
        assert!(matches!(
            "B:I,ndls".parse::<Value>(),
            Err(ParseError::InvalidIntValue(_))
        ));

        assert_eq!("B:f,0,1".parse(), Ok(Value::FloatArray(vec![0.0, 1.0])));
        assert!(matches!(
            "B:f,".parse::<Value>(),
            Err(ParseError::InvalidFloatValue(_))
        ));
        assert!(matches!(
            "B:f,ndls".parse::<Value>(),
            Err(ParseError::InvalidFloatValue(_))
        ));
    }
}
