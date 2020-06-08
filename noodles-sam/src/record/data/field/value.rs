pub mod subtype;
pub mod ty;

use std::{error, fmt, num, str::FromStr};

use super::DELIMITER;

use self::{subtype::Subtype, ty::Type};

const ARRAY_VALUE_DELIMITER: char = ',';

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Char(char),
    Int32(i32),
    Float(f32),
    String(String),
    Hex(String),
    Int8Array(Vec<i8>),
    UInt8Array(Vec<u8>),
    Int16Array(Vec<i16>),
    UInt16Array(Vec<u16>),
    Int32Array(Vec<i32>),
    UInt32Array(Vec<u32>),
    FloatArray(Vec<f32>),
}

impl Value {
    pub fn ty(&self) -> Type {
        match *self {
            Self::Char(_) => Type::Char,
            Self::Int32(_) => Type::Int32,
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

    pub fn as_char(&self) -> Option<char> {
        match *self {
            Self::Char(c) => Some(c),
            _ => None,
        }
    }

    pub fn is_char(&self) -> bool {
        self.as_char().is_some()
    }

    pub fn as_int32(&self) -> Option<i32> {
        match *self {
            Self::Int32(i) => Some(i),
            _ => None,
        }
    }

    pub fn is_int32(&self) -> bool {
        self.as_int32().is_some()
    }

    pub fn as_float(&self) -> Option<f32> {
        match *self {
            Self::Float(s) => Some(s),
            _ => None,
        }
    }

    pub fn is_float(&self) -> bool {
        self.as_float().is_some()
    }

    pub fn as_str(&self) -> Option<&str> {
        match *self {
            Self::String(ref s) => Some(s),
            _ => None,
        }
    }

    pub fn is_str(&self) -> bool {
        self.as_str().is_some()
    }

    pub fn as_hex(&self) -> Option<&str> {
        match *self {
            Self::Hex(ref h) => Some(h),
            _ => None,
        }
    }

    pub fn is_hex(&self) -> bool {
        self.as_hex().is_some()
    }

    pub fn as_int8_array(&self) -> Option<&[i8]> {
        match *self {
            Self::Int8Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_int8_array(&self) -> bool {
        self.as_int8_array().is_some()
    }

    pub fn as_uint8_array(&self) -> Option<&[u8]> {
        match *self {
            Self::UInt8Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_uint8_array(&self) -> bool {
        self.as_uint8_array().is_some()
    }

    pub fn as_int16_array(&self) -> Option<&[i16]> {
        match *self {
            Self::Int16Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_int16_array(&self) -> bool {
        self.as_int16_array().is_some()
    }

    pub fn as_uint16_array(&self) -> Option<&[u16]> {
        match *self {
            Self::UInt16Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_uint16_array(&self) -> bool {
        self.as_uint16_array().is_some()
    }

    pub fn as_int32_array(&self) -> Option<&[i32]> {
        match *self {
            Self::Int32Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_int32_array(&self) -> bool {
        self.as_int32_array().is_some()
    }

    pub fn as_uint32_array(&self) -> Option<&[u32]> {
        match *self {
            Self::UInt32Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_uint32_array(&self) -> bool {
        self.as_uint32_array().is_some()
    }

    pub fn as_float_array(&self) -> Option<&[f32]> {
        match *self {
            Self::FloatArray(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_float_array(&self) -> bool {
        self.as_float_array().is_some()
    }
}

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Char(c) => write!(f, "{}", c),
            Self::Int32(n) => write!(f, "{}", n),
            Self::Float(n) => write!(f, "{}", n),
            Self::String(s) => f.write_str(s),
            Self::Hex(s) => f.write_str(s),
            Self::Int8Array(values) => {
                f.write_str("c")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
            Self::UInt8Array(values) => {
                f.write_str("C")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
            Self::Int16Array(values) => {
                f.write_str("s")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
            Self::UInt16Array(values) => {
                f.write_str("S")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
            Self::Int32Array(values) => {
                f.write_str("i")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
            Self::UInt32Array(values) => {
                f.write_str("I")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
            Self::FloatArray(values) => {
                f.write_str("f")?;

                for value in values {
                    write!(f, ",{}", value)?
                }

                Ok(())
            }
        }
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingType,
    InvalidType(ty::ParseError),
    MissingValue,
    InvalidCharValue,
    InvalidIntValue(num::ParseIntError),
    InvalidFloatValue(num::ParseFloatError),
    MissingSubtype,
    InvalidSubtype(subtype::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingType => f.write_str("missing type"),
            Self::InvalidType(e) => write!(f, "{}", e),
            Self::MissingValue => f.write_str("missing value"),
            Self::InvalidCharValue => f.write_str("invalid char value"),
            Self::InvalidIntValue(e) => write!(f, "{}", e),
            Self::InvalidFloatValue(e) => write!(f, "{}", e),
            Self::MissingSubtype => f.write_str("missing subtype"),
            Self::InvalidSubtype(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Value {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut components = s.splitn(2, DELIMITER);

        let ty = components
            .next()
            .ok_or_else(|| ParseError::MissingType)
            .and_then(|t| t.parse().map_err(ParseError::InvalidType))?;

        let value = components.next().ok_or_else(|| ParseError::MissingValue)?;

        match ty {
            Type::Char => parse_char(value).map(Value::Char),
            Type::Int32 => parse_i32(value).map(Value::Int32),
            Type::Float => parse_f32(value).map(Value::Float),
            Type::String => Ok(Value::String(value.into())),
            Type::Hex => Ok(Value::Hex(value.into())),
            Type::Array => parse_array(value),
        }
    }
}

fn parse_char(s: &str) -> Result<char, ParseError> {
    s.chars().next().ok_or_else(|| ParseError::InvalidCharValue)
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

fn parse_f32(s: &str) -> Result<f32, ParseError> {
    s.parse().map_err(ParseError::InvalidFloatValue)
}

fn parse_array(s: &str) -> Result<Value, ParseError> {
    let mut raw_values = s.split(ARRAY_VALUE_DELIMITER);

    let subtype = raw_values
        .next()
        .ok_or_else(|| ParseError::MissingSubtype)
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
        assert_eq!(Value::Int32(0).ty(), Type::Int32);
        assert_eq!(Value::Float(0.0).ty(), Type::Float);
        assert_eq!(Value::String(String::from("noodles")).ty(), Type::String);
        assert_eq!(Value::Hex(String::from("cafe")).ty(), Type::Hex);
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
        assert_eq!(Value::Int32(0).subtype(), None);
        assert_eq!(Value::Float(0.0).subtype(), None);
        assert_eq!(Value::String(String::from("noodles")).subtype(), None);
        assert_eq!(Value::Hex(String::from("cafe")).subtype(), None);
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
        assert_eq!(Value::Int32(13).to_string(), "13");
        assert_eq!(Value::Float(3.14).to_string(), "3.14");

        assert_eq!(
            Value::String(String::from("noodles")).to_string(),
            "noodles"
        );

        assert_eq!(Value::Hex(String::from("cafe")).to_string(), "cafe");
        assert_eq!(Value::Int8Array(vec![1, -2]).to_string(), "c,1,-2");
        assert_eq!(Value::UInt8Array(vec![3, 5]).to_string(), "C,3,5");
        assert_eq!(Value::Int16Array(vec![8, -13]).to_string(), "s,8,-13");
        assert_eq!(Value::UInt16Array(vec![21, 34]).to_string(), "S,21,34");
        assert_eq!(Value::Int32Array(vec![55, -89]).to_string(), "i,55,-89");
        assert_eq!(Value::UInt32Array(vec![144, 233]).to_string(), "I,144,233");

        assert_eq!(
            Value::FloatArray(vec![2.71, 3.14]).to_string(),
            "f,2.71,3.14"
        );
    }
}
