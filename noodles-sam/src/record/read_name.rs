use std::{error, fmt, ops::Deref, str::FromStr};

use super::NULL_FIELD;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadName(Option<String>);

impl ReadName {
    pub fn is_empty(&self) -> bool {
        self.0.is_none()
    }
}

impl AsRef<str> for ReadName {
    fn as_ref(&self) -> &str {
        self.0.as_deref().unwrap_or(NULL_FIELD)
    }
}

impl Deref for ReadName {
    type Target = Option<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ReadName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid read name: {}", self.0)
    }
}

impl FromStr for ReadName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError(s.into()));
        }

        match s {
            NULL_FIELD => Ok(ReadName(None)),
            _ => Ok(ReadName(Some(s.into()))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() -> Result<(), ParseError> {
        let read_name = ReadName::default();
        assert!(read_name.is_empty());

        let read_name: ReadName = "r0".parse()?;
        assert!(!read_name.is_empty());

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let read_name = ReadName::default();
        assert_eq!(read_name.to_string(), "*");

        let read_name: ReadName = "r0".parse()?;
        assert_eq!(read_name.to_string(), "r0");

        Ok(())
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let read_name: ReadName = "*".parse()?;
        assert_eq!(*read_name, None);

        let read_name: ReadName = "r0".parse()?;
        assert_eq!(*read_name, Some(String::from("r0")));

        Ok(())
    }
}
