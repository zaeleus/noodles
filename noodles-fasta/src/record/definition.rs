//! FASTA record definition and components.

use std::fmt;

use bstr::{BStr, BString};

const PREFIX: char = '>';

/// A FASTA record definition.
///
/// A definition represents a definition line, i.e, a reference sequence name and, optionally, a
/// description.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Definition {
    name: BString,
    description: Option<BString>,
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
    pub fn new<N>(name: N, description: Option<BString>) -> Self
    where
        N: Into<BString>,
    {
        Self {
            name: name.into(),
            description,
        }
    }

    /// Returns the record name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Definition;
    /// let definition = Definition::new("sq0", None);
    /// assert_eq!(definition.name(), b"sq0");
    /// ```
    pub fn name(&self) -> &BStr {
        self.name.as_ref()
    }

    /// Returns the description if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::{BStr, BString};
    /// use noodles_fasta::record::Definition;
    ///
    /// let definition = Definition::new("sq0", None);
    /// assert_eq!(definition.description(), None);
    ///
    /// let definition = Definition::new("sq0", Some(BString::from("LN:13")));
    /// assert_eq!(definition.description(), Some(BStr::new("LN:13")));
    /// ```
    pub fn description(&self) -> Option<&BStr> {
        self.description.as_ref().map(|s| s.as_ref())
    }
}

impl fmt::Display for Definition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{PREFIX}{}", self.name)?;

        if let Some(description) = self.description() {
            write!(f, " {description}")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let definition = Definition::new("sq0", None);
        assert_eq!(definition.to_string(), ">sq0");

        let definition = Definition::new("sq0", Some(BString::from("LN:13")));
        assert_eq!(definition.to_string(), ">sq0 LN:13");
    }
}
