use std::{error, fmt, num};

use noodles_sam::header::ReferenceSequences;

// Position coordinates are 1-based.
const MIN_POSITION: i32 = 1;

static UNMAPPED_NAME: &str = "*";
static ALL_NAME: &str = ".";

/// A genomic region.
///
/// Genomic regions can either be mapped to a reference sequence, unmapped (*), or an inclusion of
/// all reads (.).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Region {
    /// A mapped region.
    Mapped {
        /// The reference sequence name.
        name: String,
        /// The start position of the region (1-based).
        start: i32,
        /// The end position of the region (1-based).
        end: i32,
    },
    /// An unmapped region.
    Unmapped,
    /// All reads.
    All,
}

impl Region {
    /// Parses a string to a region.
    ///
    /// A region string is specified as
    /// `<reference-sequence-name>[:<start-position>[-<end-position>]]`.
    ///
    /// The reference sequence name can be "*" to represent unmapped records; or ".", all records.
    /// Otherwise, the reference sequence name must exist in the reference sequence dictionary.
    ///
    /// If no start position is given, the minimum position of 1 is used. If no end position is
    /// given, the entire span of the reference sequence, i.e., its length, is used.
    pub fn from_str_reference_sequences(
        s: &str,
        reference_sequences: &ReferenceSequences,
    ) -> Result<Self, ParseError> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        } else if s == UNMAPPED_NAME {
            return Ok(Self::Unmapped);
        } else if s == ALL_NAME {
            return Ok(Self::All);
        }

        if let Some(i) = s.rfind(':') {
            let suffix = &s[i + 1..];

            if let Ok((start, end)) = parse_interval(suffix) {
                let prefix = &s[0..i];

                if let Some(reference_sequence) = reference_sequences.get(prefix) {
                    if reference_sequences.contains_key(s) {
                        return Err(ParseError::Ambiguous);
                    } else {
                        let resolved_end = end.unwrap_or(reference_sequence.len() as i32);
                        return Ok(Self::mapped(prefix, start, resolved_end));
                    }
                }
            }
        }

        if let Some(reference_sequence) = reference_sequences.get(s) {
            let end = reference_sequence.len() as i32;
            Ok(Self::mapped(s, MIN_POSITION, end))
        } else {
            Err(ParseError::Invalid)
        }
    }

    /// Creates a new mapped region.
    ///
    /// `start` and `end` are the start and (optional) end positions of the region in the given
    /// reference sequence `name`. When `end` is `None`, most analyses will assume the end to be
    /// unbounded, i.e., until the end of the reference sequence.
    ///
    /// Positions are assumed to be 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    /// let region = Region::mapped("sq0", 1, 5);
    /// assert!(matches!(region, Region::Mapped { name, start: 1, end: 5 }));
    /// ```
    pub fn mapped<I>(name: I, start: i32, end: i32) -> Self
    where
        I: Into<String>,
    {
        Self::Mapped {
            name: name.into(),
            start,
            end,
        }
    }

    /// Returns the reference name of the region.
    ///
    /// If the region is unmapped, this returns "*". If the region represents
    /// all records, this returns ".".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Region;
    ///
    /// let region = Region::mapped("sq0", 1, 5);
    /// assert_eq!(region.name(), "sq0");
    ///
    /// assert_eq!(Region::Unmapped.name(), "*");
    /// assert_eq!(Region::All.name(), ".");
    /// ```
    pub fn name(&self) -> &str {
        match self {
            Self::Mapped { name, .. } => name,
            Self::Unmapped => UNMAPPED_NAME,
            Self::All => ALL_NAME,
        }
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Mapped { name, start, end } => write!(f, "{}:{}-{}", name, start, end),
            Self::Unmapped => write!(f, "{}", UNMAPPED_NAME),
            Self::All => write!(f, "{}", ALL_NAME),
        }
    }
}

/// An error returned when a genomic region fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is ambiguous.
    Ambiguous,
    /// The input is invalid.
    Invalid,
    /// The start position is invalid.
    InvalidStartPosition(num::ParseIntError),
    /// The end position is invalid.
    InvalidEndPosition(num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Ambiguous => f.write_str("ambiguous input"),
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidStartPosition(e) => write!(f, "invalid start position: {}", e),
            Self::InvalidEndPosition(e) => write!(f, "invalid end position: {}", e),
        }
    }
}

fn parse_interval(s: &str) -> Result<(i32, Option<i32>), ParseError> {
    let mut components = s.splitn(2, '-');

    let start = match components.next() {
        Some(t) => t.parse().map_err(ParseError::InvalidStartPosition)?,
        None => MIN_POSITION,
    };

    let end = match components.next() {
        Some(t) => t
            .parse()
            .map(Some)
            .map_err(ParseError::InvalidEndPosition)?,
        None => None,
    };

    Ok((start, end))
}

#[cfg(test)]
mod tests {
    use noodles_sam::header::ReferenceSequence;

    use super::*;

    #[test]
    fn test_from_str_reference_sequences() {
        let reference_sequences: ReferenceSequences = vec![
            ReferenceSequence::new(String::from("sq0"), 8),
            ReferenceSequence::new(String::from("sq1:"), 13),
            ReferenceSequence::new(String::from("sq2:5"), 21),
            ReferenceSequence::new(String::from("sq3"), 34),
            ReferenceSequence::new(String::from("sq3:5-8"), 55),
        ]
        .into_iter()
        .map(|rs| (rs.name().into(), rs))
        .collect();

        assert_eq!(
            Region::from_str_reference_sequences("*", &reference_sequences),
            Ok(Region::Unmapped)
        );

        assert_eq!(
            Region::from_str_reference_sequences(".", &reference_sequences),
            Ok(Region::All)
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq0:3-5", &reference_sequences),
            Ok(Region::Mapped {
                name: String::from("sq0"),
                start: 3,
                end: 5
            })
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq0:3", &reference_sequences),
            Ok(Region::Mapped {
                name: String::from("sq0"),
                start: 3,
                end: 8
            })
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq0", &reference_sequences),
            Ok(Region::Mapped {
                name: String::from("sq0"),
                start: 1,
                end: 8
            })
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq1:", &reference_sequences),
            Ok(Region::Mapped {
                name: String::from("sq1:"),
                start: 1,
                end: 13
            })
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq2:5", &reference_sequences),
            Ok(Region::Mapped {
                name: String::from("sq2:5"),
                start: 1,
                end: 21
            })
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq3:8-13", &reference_sequences),
            Ok(Region::Mapped {
                name: String::from("sq3"),
                start: 8,
                end: 13
            })
        );

        assert_eq!(
            Region::from_str_reference_sequences("sq3:5-8", &reference_sequences),
            Err(ParseError::Ambiguous)
        );

        assert_eq!(
            Region::from_str_reference_sequences("", &reference_sequences),
            Err(ParseError::Empty)
        );
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Region::mapped("sq0", 3, 5).to_string(), "sq0:3-5");
        assert_eq!(Region::Unmapped.to_string(), "*");
        assert_eq!(Region::All.to_string(), ".");
    }
}
