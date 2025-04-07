//! GFF directive genome build.

use std::{error, fmt, str::FromStr};

use bstr::{BStr, BString};

/// A GFF directive genome build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GenomeBuild {
    source: BString,
    name: BString,
}

impl GenomeBuild {
    /// Creates a genome build directive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::GenomeBuild;
    /// let genome_build = GenomeBuild::new("NDLS", "r1");
    /// ```
    pub fn new<S, N>(source: S, name: N) -> Self
    where
        S: Into<BString>,
        N: Into<BString>,
    {
        Self {
            source: source.into(),
            name: name.into(),
        }
    }

    /// Returns the genome build source.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::GenomeBuild;
    /// let genome_build = GenomeBuild::new("NDLS", "r1");
    /// assert_eq!(genome_build.source(), "NDLS");
    /// ```
    pub fn source(&self) -> &BStr {
        self.source.as_ref()
    }

    /// Returns the genome build name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::GenomeBuild;
    /// let genome_build = GenomeBuild::new("NDLS", "r1");
    /// assert_eq!(genome_build.name(), "r1");
    /// ```
    pub fn name(&self) -> &BStr {
        self.name.as_ref()
    }
}

impl fmt::Display for GenomeBuild {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.source, self.name)
    }
}

/// An error returned when a raw GFF genome build directive fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The genome build source is missing.
    MissingSource,
    /// The genome build name is missing.
    MissingName,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid genome build directive: ")?;

        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingSource => f.write_str("missing source"),
            Self::MissingName => f.write_str("missing name"),
        }
    }
}

impl FromStr for GenomeBuild {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut args = s.split_ascii_whitespace();
        let source = args.next().ok_or(ParseError::MissingSource)?;
        let name = args.next().ok_or(ParseError::MissingName)?;

        Ok(Self::new(source, name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let genome_build = GenomeBuild::new("NDLS", "r1");
        assert_eq!(genome_build.to_string(), "NDLS r1");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("NDLS r1".parse(), Ok(GenomeBuild::new("NDLS", "r1")));

        assert_eq!("".parse::<GenomeBuild>(), Err(ParseError::Empty));
        assert_eq!("NDLS".parse::<GenomeBuild>(), Err(ParseError::MissingName));

        Ok(())
    }
}
