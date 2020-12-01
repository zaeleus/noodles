//! FASTA record definition and components.

use std::{error, fmt, str::FromStr};

const PREFIX: char = '>';

/// A FASTA record definition.
///
/// A definition represents a definition line, i.e, a reference sequence name and, optionally, a
/// description.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Definition {
    reference_sequence_name: String,
    description: Option<String>,
}

impl Definition {
    /// Creates a FASTA record definition.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    /// let definition = Definition::new(String::from("sq0"), None);
    /// ```
    pub fn new(reference_sequence_name: String, description: Option<String>) -> Self {
        Self {
            reference_sequence_name,
            description,
        }
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    /// let definition = Definition::new(String::from("sq0"), None);
    /// assert_eq!(definition.reference_sequence_name(), "sq0");
    /// ```
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    /// Returns the description if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    ///
    /// let definition = Definition::new(String::from("sq0"), None);
    /// assert_eq!(definition.description(), None);
    ///
    /// let definition = Definition::new(String::from("sq0"), Some(String::from("LN:13")));
    /// assert_eq!(definition.description(), Some("LN:13"));
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
}

/// An error returned when a raw record definition fails to parse.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The prefix (`>`) is missing.
    MissingPrefix,
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("record header is empty"),
            Self::MissingPrefix => f.write_str("the '>' prefix is missing"),
            Self::MissingReferenceSequenceName => f.write_str("the reference sequence is missing"),
        }
    }
}

impl FromStr for Definition {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        } else if !s.starts_with(PREFIX) {
            return Err(ParseError::MissingPrefix);
        }

        let line = &s[1..];
        let mut components = line.splitn(2, |c: char| c.is_ascii_whitespace());

        let reference_sequence_name = components
            .next()
            .and_then(|s| if s.is_empty() { None } else { Some(s.into()) })
            .ok_or(ParseError::MissingReferenceSequenceName)?;

        let description = components.next().map(|s| s.trim().into());

        Ok(Self {
            reference_sequence_name,
            description,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let definition: Definition = ">sq0".parse()?;
        assert_eq!(definition.reference_sequence_name(), "sq0");
        assert!(definition.description().is_none());

        let definition: Definition = ">sq0  LN:13".parse()?;
        assert_eq!(definition.reference_sequence_name(), "sq0");
        assert_eq!(definition.description(), Some("LN:13"));

        assert_eq!("".parse::<Definition>(), Err(ParseError::Empty));
        assert_eq!("sq0".parse::<Definition>(), Err(ParseError::MissingPrefix));
        assert_eq!(
            ">".parse::<Definition>(),
            Err(ParseError::MissingReferenceSequenceName)
        );

        Ok(())
    }
}
