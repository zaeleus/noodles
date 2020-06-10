use std::{error, fmt, ops::Deref, str::FromStr};

use super::NULL_FIELD;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReferenceSequenceName(Option<String>);

impl ReferenceSequenceName {
    pub fn is_empty(&self) -> bool {
        self.is_none()
    }
}

impl AsRef<str> for ReferenceSequenceName {
    fn as_ref(&self) -> &str {
        self.0.as_deref().unwrap_or(NULL_FIELD)
    }
}

impl Deref for ReferenceSequenceName {
    type Target = Option<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ReferenceSequenceName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

#[derive(Debug)]
pub enum ParseError {
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("reference sequence name cannot be empty"),
        }
    }
}

impl FromStr for ReferenceSequenceName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            NULL_FIELD => Ok(ReferenceSequenceName(None)),
            _ => Ok(ReferenceSequenceName(Some(s.into()))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let reference_sequence_name = ReferenceSequenceName::default();
        assert!(reference_sequence_name.is_empty());

        let reference_sequence_name = ReferenceSequenceName(Some(String::from("sq0")));
        assert!(!reference_sequence_name.is_empty());
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let reference_sequence_name = ReferenceSequenceName::default();
        assert_eq!(reference_sequence_name.to_string(), "*");

        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        assert_eq!(reference_sequence_name.to_string(), "sq0");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let reference_sequence_name: ReferenceSequenceName = "*".parse()?;
        assert_eq!(*reference_sequence_name, None);

        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        assert_eq!(*reference_sequence_name, Some(String::from("sq0")));

        assert!("".parse::<ReferenceSequenceName>().is_err());

        Ok(())
    }
}
