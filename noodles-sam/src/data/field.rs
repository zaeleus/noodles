mod value;

pub use value::Value;

use std::{error, fmt, str::FromStr};

pub(crate) const DELIMITER: char = ':';
const TAG_LEN: usize = 2;

#[derive(Clone, Debug)]
pub struct Field {
    tag: String,
    value: Value,
}

impl Field {
    pub fn new(tag: String, value: Value) -> Self {
        Self { tag, value }
    }

    pub fn tag(&self) -> &str {
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

#[derive(Clone, Copy, Debug)]
pub enum Component {
    Tag,
    Type,
    Subtype,
    Value,
}

#[derive(Debug)]
pub enum ParseError {
    Missing(Component),
    Invalid(Component, Box<dyn std::error::Error>),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(component) => write!(f, "missing field component: {:?}", component),
            Self::Invalid(component, e) => {
                write!(f, "invalid field {:?} component: {}", component, e)
            }
        }
    }
}

impl FromStr for Field {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut components = s.splitn(2, DELIMITER);

        let tag = components
            .next()
            .ok_or_else(|| ParseError::Missing(Component::Tag))
            .and_then(parse_tag)?;

        let value = components
            .next()
            .ok_or_else(|| ParseError::Missing(Component::Type))
            .and_then(|t| t.parse())?;

        Ok(Self::new(tag.into(), value))
    }
}

fn parse_tag(s: &str) -> Result<&str, ParseError> {
    if s.len() == TAG_LEN {
        Ok(s)
    } else {
        Err(ParseError::Invalid(
            Component::Tag,
            format!("invalid tag length: expected {}, got {}", TAG_LEN, s.len()).into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let field = Field::new(String::from("RG"), Value::String(String::from("rg0")));
        assert_eq!(field.to_string(), "RG:Z:rg0");
    }

    #[test]
    fn test_parse_tag() {
        assert_eq!(parse_tag("RG").unwrap(), "RG");
        assert!(parse_tag("").is_err());
        assert!(parse_tag("R").is_err());
        assert!(parse_tag("noodles").is_err());
    }
}
