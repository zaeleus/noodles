//! GFF directive sequence region.

use std::{error, fmt, num, str::FromStr};

use bstr::{BStr, BString};

/// A GFF directive sequence region.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceRegion {
    reference_sequence_name: BString,
    start: i32,
    end: i32,
}

impl SequenceRegion {
    /// Creates a sequence region directive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::SequenceRegion;
    /// let sequence_region = SequenceRegion::new("sq0", 8, 13);
    /// ```
    pub fn new<N>(reference_sequence_name: N, start: i32, end: i32) -> Self
    where
        N: Into<BString>,
    {
        Self {
            reference_sequence_name: reference_sequence_name.into(),
            start,
            end,
        }
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::SequenceRegion;
    /// let sequence_region = SequenceRegion::new("sq0", 8, 13);
    /// assert_eq!(sequence_region.reference_sequence_name(), "sq0");
    /// ```
    pub fn reference_sequence_name(&self) -> &BStr {
        self.reference_sequence_name.as_ref()
    }

    /// Returns the start position.
    ///
    /// The start position is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::SequenceRegion;
    /// let sequence_region = SequenceRegion::new("sq0", 8, 13);
    /// assert_eq!(sequence_region.start(), 8);
    /// ```
    pub fn start(&self) -> i32 {
        self.start
    }

    /// Returns the end position.
    ///
    /// The end position is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff::directive_buf::value::SequenceRegion;
    /// let sequence_region = SequenceRegion::new("sq0", 8, 13);
    /// assert_eq!(sequence_region.end(), 13);
    /// ```
    pub fn end(&self) -> i32 {
        self.end
    }
}

impl fmt::Display for SequenceRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} {}",
            self.reference_sequence_name, self.start, self.end
        )
    }
}

/// An error returned when a raw GFF sequence region directive fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The start is missing.
    MissingStart,
    /// The start is invalid.
    InvalidStart(num::ParseIntError),
    /// The end is missing.
    MissingEnd,
    /// The end is invalid.
    InvalidEnd(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStart(e) | Self::InvalidEnd(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingReferenceSequenceName => f.write_str("missing reference sequence name"),
            Self::MissingStart => f.write_str("missing start"),
            Self::InvalidStart(_) => f.write_str("invalid start"),
            Self::MissingEnd => f.write_str("missing end"),
            Self::InvalidEnd(_) => f.write_str("invalid end"),
        }
    }
}

impl FromStr for SequenceRegion {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut args = s.split_ascii_whitespace();

        let reference_sequence_name = args
            .next()
            .ok_or(ParseError::MissingReferenceSequenceName)?;

        let start = args
            .next()
            .ok_or(ParseError::MissingStart)
            .and_then(|s| s.parse().map_err(ParseError::InvalidStart))?;

        let end = args
            .next()
            .ok_or(ParseError::MissingEnd)
            .and_then(|s| s.parse().map_err(ParseError::InvalidEnd))?;

        Ok(Self::new(reference_sequence_name, start, end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let sequence_region = SequenceRegion::new("sq0", 8, 13);
        assert_eq!(sequence_region.to_string(), "sq0 8 13");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!(
            "sq0 8 13".parse::<SequenceRegion>()?,
            SequenceRegion::new("sq0", 8, 13)
        );

        assert_eq!("".parse::<SequenceRegion>(), Err(ParseError::Empty));

        assert_eq!(
            "sq0".parse::<SequenceRegion>(),
            Err(ParseError::MissingStart)
        );

        assert!(matches!(
            "sq0 eight".parse::<SequenceRegion>(),
            Err(ParseError::InvalidStart(_))
        ));

        assert_eq!(
            "sq0 8".parse::<SequenceRegion>(),
            Err(ParseError::MissingEnd)
        );

        assert!(matches!(
            "sq0 8 thirteen".parse::<SequenceRegion>(),
            Err(ParseError::InvalidEnd(_))
        ));

        Ok(())
    }
}
