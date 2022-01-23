//! FASTA record definition and components.

use std::{error, fmt, str::FromStr};

const PREFIX: char = '>';

/// A FASTA record definition.
///
/// A definition represents a definition line, i.e, a reference sequence name and, optionally, a
/// description.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Definition {
    name: String,
    description: Option<String>,
}

impl Definition {
    /// Creates a FASTA record definition.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    /// let definition = Definition::new("sq0", None);
    /// ```
    pub fn new<N>(name: N, description: Option<String>) -> Self
    where
        N: Into<String>,
    {
        Self {
            name: name.into(),
            description,
        }
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    /// let definition = Definition::new("sq0", None);
    /// assert_eq!(definition.reference_sequence_name(), "sq0");
    /// ```
    #[deprecated(since = "0.3.0", note = "Use `name` instead.")]
    pub fn reference_sequence_name(&self) -> &str {
        &self.name
    }

    /// Returns the record name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    /// let definition = Definition::new("sq0", None);
    /// assert_eq!(definition.name(), "sq0");
    /// ```
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the description if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    ///
    /// let definition = Definition::new("sq0", None);
    /// assert_eq!(definition.description(), None);
    ///
    /// let definition = Definition::new("sq0", Some(String::from("LN:13")));
    /// assert_eq!(definition.description(), Some("LN:13"));
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
}

impl fmt::Display for Definition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", PREFIX, self.name())?;

        if let Some(description) = self.description() {
            write!(f, " {}", description)?;
        }

        Ok(())
    }
}

/// An error returned when a raw record definition fails to parse.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The prefix (`>`) is missing.
    MissingPrefix,
    /// The sequence name is missing.
    MissingName,
    /// The sequence name is missing.
    #[deprecated(since = "0.6.0", note = "Use `ParseError::MissingName` instead.")]
    MissingReferenceSequenceName,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingPrefix => write!(f, "missing prefix ('{}')", PREFIX),
            Self::MissingName => f.write_str("missing name"),
            #[allow(deprecated)]
            Self::MissingReferenceSequenceName => f.write_str("missing name"),
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

        let name = components
            .next()
            .and_then(|s| if s.is_empty() { None } else { Some(s.into()) })
            .ok_or(ParseError::MissingName)?;

        let description = components.next().map(|s| s.trim().into());

        Ok(Self { name, description })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let definition = Definition::new("sq0", None);
        assert_eq!(definition.to_string(), ">sq0");

        let definition = Definition::new("sq0", Some(String::from("LN:13")));
        assert_eq!(definition.to_string(), ">sq0 LN:13");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(">sq0".parse(), Ok(Definition::new("sq0", None)));

        assert_eq!(
            ">sq0  LN:13".parse(),
            Ok(Definition::new("sq0", Some(String::from("LN:13"))))
        );

        assert_eq!("".parse::<Definition>(), Err(ParseError::Empty));
        assert_eq!("sq0".parse::<Definition>(), Err(ParseError::MissingPrefix));
        assert_eq!(">".parse::<Definition>(), Err(ParseError::MissingName));
    }
}
