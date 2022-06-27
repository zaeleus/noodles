//! SAM record data field and components.

pub mod tag;
pub mod value;

pub use self::{tag::Tag, value::Value};

use std::{error, fmt, str::FromStr};

pub(crate) const DELIMITER: char = ':';

/// A SAM record data field.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    tag: Tag,
    value: Value,
}

impl Field {
    /// Creates a SAM record data field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::{field::{Tag, Value}, Field};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// ```
    pub fn new(tag: Tag, value: Value) -> Self {
        Self { tag, value }
    }

    /// Returns the data field tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::{field::{Tag, Value}, Field};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// assert_eq!(field.tag(), Tag::AlignmentHitCount);
    /// ```
    pub fn tag(&self) -> Tag {
        self.tag
    }

    /// Returns the data field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::data::{field::{Tag, Value}, Field};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// assert!(matches!(field.value(), Value::Int32(n) if *n == 1));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use self::value::Type;

        let ty = if self.value().is_int() {
            Type::Int32
        } else {
            self.value().ty()
        };

        write!(
            f,
            "{}{}{}{}{}",
            self.tag, DELIMITER, ty, DELIMITER, self.value
        )
    }
}

/// An error returned when a raw SAM record data field fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The data field tag is invalid.
    InvalidTag(tag::ParseError),
    /// The data field type is invalid.
    InvalidType(value::ty::ParseError),
    /// The data field value is invalid.
    InvalidValue,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidTag(e) => write!(f, "invalid tag: {}", e),
            Self::InvalidType(e) => write!(f, "invalid type: {}", e),
            Self::InvalidValue => f.write_str("invalid value"),
        }
    }
}

impl FromStr for Field {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (t, rest) = s.split_once(DELIMITER).ok_or(ParseError::Invalid)?;
        let tag = t.parse().map_err(ParseError::InvalidTag)?;

        let (raw_ty, raw_value) = rest.split_once(DELIMITER).ok_or(ParseError::Invalid)?;
        let ty = raw_ty.parse().map_err(ParseError::InvalidType)?;
        let value = Value::from_str_type(raw_value, ty).map_err(|_| ParseError::InvalidValue)?;

        Ok(Self::new(tag, value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let field = Field::new(Tag::AlignmentHitCount, Value::UInt8(1));
        assert_eq!(field.to_string(), "NH:i:1");

        let field = Field::new(Tag::ReadGroup, Value::String(String::from("rg0")));
        assert_eq!(field.to_string(), "RG:Z:rg0");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "RG:Z:rg0".parse(),
            Ok(Field::new(
                Tag::ReadGroup,
                Value::String(String::from("rg0"))
            ))
        );

        assert_eq!("".parse::<Field>(), Err(ParseError::Invalid));

        assert!(matches!(
            "_:Z:rg0".parse::<Field>(),
            Err(ParseError::InvalidTag(_))
        ));

        assert!(matches!(
            "RG:_:rg0".parse::<Field>(),
            Err(ParseError::InvalidType(_))
        ));

        assert_eq!("RG:Z:🍜".parse::<Field>(), Err(ParseError::InvalidValue));
    }
}
