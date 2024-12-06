//! GFF directives.

pub mod key;
pub mod value;

use std::{error, fmt, str::FromStr};

use crate::Directive;

pub use self::value::Value;

pub(crate) const PREFIX: &str = "##";

/// A GFF directive.
///
/// This is also called a pragma or metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DirectiveBuf {
    key: String,
    value: Option<Value>,
}

impl DirectiveBuf {
    /// Creates a directive buffer.
    pub fn new<K>(key: K, value: Option<Value>) -> Self
    where
        K: Into<String>,
    {
        Self {
            key: key.into(),
            value,
        }
    }

    /// Returns the key.
    pub fn key(&self) -> &str {
        &self.key
    }

    /// Returns the value.
    pub fn value(&self) -> Option<&Value> {
        self.value.as_ref()
    }
}

impl From<Directive<'_>> for DirectiveBuf {
    fn from(directive: Directive<'_>) -> Self {
        Self::new(directive.key(), directive.value().map(Value::from))
    }
}

/// An error returned when a raw GFF directive fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The directive prefix (`##`) is missing.
    MissingPrefix,
    /// The directive name is missing.
    MissingKey,
    /// The directive value is missing.
    MissingValue,
    /// The GFF version is invalid.
    InvalidGffVersion(value::gff_version::ParseError),
    /// A sequence region is invalid.
    InvalidSequenceRegion(value::sequence_region::ParseError),
    /// A genome build is invalid.
    InvalidGenomeBuild(value::genome_build::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidGffVersion(e) => Some(e),
            Self::InvalidSequenceRegion(e) => Some(e),
            Self::InvalidGenomeBuild(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPrefix => f.write_str("directive prefix is missing"),
            Self::MissingKey => f.write_str("directive key is missing"),
            Self::MissingValue => f.write_str("directive value is missing"),
            Self::InvalidGffVersion(_) => f.write_str("invalid GFF version"),
            Self::InvalidSequenceRegion(_) => f.write_str("invalid sequence region"),
            Self::InvalidGenomeBuild(_) => f.write_str("invalid genome build"),
        }
    }
}

impl FromStr for DirectiveBuf {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !s.starts_with(PREFIX) {
            return Err(ParseError::MissingPrefix);
        }

        let mut components = s[PREFIX.len()..].splitn(2, |c: char| c.is_ascii_whitespace());

        let key = components.next().ok_or(ParseError::MissingKey)?;

        let value = match key {
            key::GFF_VERSION => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidGffVersion))
                .map(Value::GffVersion)
                .map(Some)?,
            key::SEQUENCE_REGION => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidSequenceRegion))
                .map(Value::SequenceRegion)
                .map(Some)?,
            key::GENOME_BUILD => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidGenomeBuild))
                .map(Value::GenomeBuild)
                .map(Some)?,
            _ => components.next().map(Value::from),
        };

        Ok(Self {
            key: key.into(),
            value,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("##noodles".parse(), Ok(DirectiveBuf::new("noodles", None)));

        assert_eq!(
            "##noodles gff".parse(),
            Ok(DirectiveBuf::new("noodles", Some(Value::from("gff")))),
        );
    }
}
