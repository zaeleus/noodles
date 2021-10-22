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
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int(1));
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
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int(1));
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
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int(1));
    /// assert!(matches!(field.value(), Value::Int(n) if *n == 1));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}{}{}{}",
            self.tag,
            DELIMITER,
            self.value.ty(),
            DELIMITER,
            self.value
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
    /// The data field value is invalid.
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidTag(e) => write!(f, "invalid tag: {}", e),
            Self::InvalidValue(e) => write!(f, "invalid value: {}", e),
        }
    }
}

impl FromStr for Field {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.split_once(DELIMITER) {
            Some((t, v)) => {
                let tag = t.parse().map_err(ParseError::InvalidTag)?;
                let value = v.parse().map_err(ParseError::InvalidValue)?;
                Ok(Self::new(tag, value))
            }
            None => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let field = Field::new(Tag::ReadGroup, Value::String(String::from("rg0")));
        assert_eq!(field.to_string(), "RG:Z:rg0");
    }
}
