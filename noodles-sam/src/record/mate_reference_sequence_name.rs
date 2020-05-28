use std::{error, fmt, str::FromStr};

use super::NULL_FIELD;

const EQ_FIELD: &str = "=";

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MateReferenceSequenceName {
    None,
    Eq,
    Some(String),
}

impl MateReferenceSequenceName {
    pub fn as_str(&self) -> &str {
        match self {
            Self::None => NULL_FIELD,
            Self::Eq => EQ_FIELD,
            Self::Some(name) => name,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.is_none()
    }

    pub fn is_none(&self) -> bool {
        match self {
            Self::None => true,
            _ => false,
        }
    }

    pub fn is_eq(&self) -> bool {
        match self {
            Self::Eq => true,
            _ => false,
        }
    }

    pub fn is_some(&self) -> bool {
        match self {
            Self::Some(_) => true,
            _ => false,
        }
    }
}

impl Default for MateReferenceSequenceName {
    fn default() -> Self {
        Self::None
    }
}

impl fmt::Display for MateReferenceSequenceName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid mate reference sequence name: {}", self.0)
    }
}

impl FromStr for MateReferenceSequenceName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            EQ_FIELD => Ok(MateReferenceSequenceName::Eq),
            NULL_FIELD => Ok(MateReferenceSequenceName::None),
            _ => Ok(MateReferenceSequenceName::Some(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_str() {
        assert_eq!(MateReferenceSequenceName::None.as_str(), "*");
        assert_eq!(MateReferenceSequenceName::Eq.as_str(), "=");
        assert_eq!(
            MateReferenceSequenceName::Some(String::from("sq0")).as_str(),
            "sq0"
        )
    }

    #[test]
    fn test_is_empty() {
        assert!(MateReferenceSequenceName::None.is_empty());
        assert!(!MateReferenceSequenceName::Eq.is_empty());
        assert!(!MateReferenceSequenceName::Some(String::from("sq0")).is_empty());
    }

    #[test]
    fn test_is_none() {
        assert!(MateReferenceSequenceName::None.is_none());
        assert!(!MateReferenceSequenceName::Eq.is_none());
        assert!(!MateReferenceSequenceName::Some(String::from("sq0")).is_none());
    }

    #[test]
    fn test_is_eq() {
        assert!(!MateReferenceSequenceName::None.is_eq());
        assert!(MateReferenceSequenceName::Eq.is_eq());
        assert!(!MateReferenceSequenceName::Some(String::from("sq0")).is_eq());
    }

    #[test]
    fn test_is_some() {
        assert!(!MateReferenceSequenceName::None.is_some());
        assert!(!MateReferenceSequenceName::Eq.is_some());
        assert!(MateReferenceSequenceName::Some(String::from("sq0")).is_some());
    }

    #[test]
    fn test_default() {
        assert!(MateReferenceSequenceName::default().is_empty());
    }
}
