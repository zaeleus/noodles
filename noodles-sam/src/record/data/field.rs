pub mod tag;
pub mod value;

pub use self::{tag::Tag, value::Value};

use std::{error, fmt, str::FromStr};

pub(crate) const DELIMITER: char = ':';

#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    tag: Tag,
    value: Value,
}

impl Field {
    pub fn new(tag: Tag, value: Value) -> Self {
        Self { tag, value }
    }

    pub fn tag(&self) -> &Tag {
        &self.tag
    }

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

#[derive(Debug)]
pub enum ParseError {
    MissingTag,
    InvalidTag(tag::ParseError),
    MissingValue,
    InvalidValue(value::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTag => f.write_str("missing tag"),
            Self::InvalidTag(e) => write!(f, "{}", e),
            Self::MissingValue => f.write_str("missing value"),
            Self::InvalidValue(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Field {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut components = s.splitn(2, DELIMITER);

        let tag = components
            .next()
            .ok_or_else(|| ParseError::MissingTag)
            .and_then(|t| t.parse().map_err(ParseError::InvalidTag))?;

        let value = components
            .next()
            .ok_or_else(|| ParseError::MissingValue)
            .and_then(|t| t.parse().map_err(ParseError::InvalidValue))?;

        Ok(Self::new(tag, value))
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
