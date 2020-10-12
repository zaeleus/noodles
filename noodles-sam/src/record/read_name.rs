//! SAM record read name.

use std::{error, fmt, ops::Deref, str::FromStr};

use super::NULL_FIELD;

// § 1.4 The alignment section: mandatory fields (2020-07-19)
const MAX_LENGTH: usize = 254;

/// A SAM record read name.
///
/// This is also called a query name.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadName(Option<String>);

impl ReadName {
    /// Returns whether a read name is set.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam as sam;
    /// use noodles_sam::record::ReadName;
    ///
    /// let read_name = ReadName::default();
    /// assert!(read_name.is_empty());
    ///
    /// let read_name: ReadName = "r0".parse()?;
    /// assert!(!read_name.is_empty());
    /// # Ok::<(), sam::record::read_name::ParseError>(())
    /// ```
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

/// An error returned when a raw SAM record read name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for ReadName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            NULL_FIELD => Ok(ReadName(None)),
            _ => {
                if s.len() > MAX_LENGTH {
                    Err(ParseError::Invalid)
                } else {
                    Ok(ReadName(Some(s.into())))
                }
            }
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
    fn test_from_str() {
        assert_eq!("*".parse(), Ok(ReadName(None)));
        assert_eq!("r0".parse(), Ok(ReadName(Some(String::from("r0")))));

        assert_eq!("".parse::<ReadName>(), Err(ParseError::Empty));

        let s: String = (0..MAX_LENGTH + 1).map(|_| 'N').collect();
        assert_eq!(s.parse::<ReadName>(), Err(ParseError::Invalid));
    }
}
