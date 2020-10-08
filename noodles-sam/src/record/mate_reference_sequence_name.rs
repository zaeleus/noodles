//! SAM record mate reference sequence name.

use std::{error, fmt, str::FromStr};

use super::NULL_FIELD;

const EQ_FIELD: &str = "=";

/// A SAM record mate reference sequence name.
///
/// SAM record mate reference sequence names can be in 1 of 3 states:
///
///   1. the mate is not set (`*`),
///   2. the mate is on the same reference sequence (`=`), or
///   3. the mate is on a different reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum MateReferenceSequenceName {
    /// The mate is not set.
    None,
    /// The mate is on the same reference sequence.
    Eq,
    /// The mate is on a different reference sequence.
    Some(String),
}

impl MateReferenceSequenceName {
    /// Returns whether the mate is not set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::MateReferenceSequenceName;
    /// assert!(MateReferenceSequenceName::None.is_empty());
    /// assert!(!MateReferenceSequenceName::Eq.is_empty());
    /// assert!(!MateReferenceSequenceName::Some(String::from("sq0")).is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.is_none()
    }

    /// Returns whether the mate is not set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::MateReferenceSequenceName;
    /// assert!(MateReferenceSequenceName::None.is_none());
    /// assert!(!MateReferenceSequenceName::Eq.is_none());
    /// assert!(!MateReferenceSequenceName::Some(String::from("sq0")).is_none());
    /// ```
    pub fn is_none(&self) -> bool {
        matches!(self, Self::None)
    }

    /// Returns whether the mate is on the same reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::MateReferenceSequenceName;
    /// assert!(!MateReferenceSequenceName::None.is_eq());
    /// assert!(MateReferenceSequenceName::Eq.is_eq());
    /// assert!(!MateReferenceSequenceName::Some(String::from("sq0")).is_eq());
    /// ```
    pub fn is_eq(&self) -> bool {
        matches!(self, Self::Eq)
    }

    /// Returns whether the mate is on a different reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::MateReferenceSequenceName;
    /// assert!(!MateReferenceSequenceName::None.is_some());
    /// assert!(!MateReferenceSequenceName::Eq.is_some());
    /// assert!(MateReferenceSequenceName::Some(String::from("sq0")).is_some());
    /// ```
    pub fn is_some(&self) -> bool {
        matches!(self, Self::Some(_))
    }
}

impl AsRef<str> for MateReferenceSequenceName {
    fn as_ref(&self) -> &str {
        match self {
            Self::None => NULL_FIELD,
            Self::Eq => EQ_FIELD,
            Self::Some(name) => name,
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
        write!(f, "{}", self.as_ref())
    }
}

/// An error returned when a raw SAM record mate reference sequence name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("mate reference sequence name cannot be empty"),
        }
    }
}

impl FromStr for MateReferenceSequenceName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
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
    fn test_default() {
        assert!(MateReferenceSequenceName::default().is_empty());
    }

    #[test]
    fn test_fmt() {
        assert_eq!(MateReferenceSequenceName::None.to_string(), "*");
        assert_eq!(MateReferenceSequenceName::Eq.to_string(), "=");
        assert_eq!(
            MateReferenceSequenceName::Some(String::from("sq0")).to_string(),
            "sq0"
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!("*".parse(), Ok(MateReferenceSequenceName::None));
        assert_eq!("=".parse(), Ok(MateReferenceSequenceName::Eq));
        assert_eq!(
            "sq0".parse(),
            Ok(MateReferenceSequenceName::Some(String::from("sq0")))
        );

        assert_eq!(
            "".parse::<MateReferenceSequenceName>(),
            Err(ParseError::Empty)
        );
    }
}
