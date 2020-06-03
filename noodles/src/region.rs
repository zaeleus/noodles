use std::{error, fmt, str::FromStr};

use noodles_sam::header::{ReferenceSequence, ReferenceSequences};

// Position coordinates are 1-based.
const MIN_POSITION: u64 = 1;

static UNMAPPED_NAME: &str = "*";
static ALL_NAME: &str = ".";

/// A genomic region.
///
/// Genomic regions can either be mapped to a reference sequence, unmapped (*), or an inclusion of
/// all reads (.).
#[derive(Clone, Debug)]
pub enum Region {
    Mapped {
        name: String,
        start: u64,
        end: Option<u64>,
    },
    Unmapped,
    All,
}

impl Region {
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
    /// use noodles::Region;
    ///
    /// let region = Region::mapped("sq0", 1, Some(5));
    ///
    /// assert!(matches!(region, Region::Mapped {
    ///     name,
    ///     start: 1,
    ///     end: Some(5),
    /// }));
    /// ```
    pub fn mapped<I>(name: I, start: u64, end: Option<u64>) -> Region
    where
        I: Into<String>,
    {
        Region::Mapped {
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
    /// use noodles::Region;
    ///
    /// let region = Region::mapped("sq0", 1, Some(5));
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

    /// Resolves the region by finding it in a given list of reference sequences.
    ///
    /// If the region name exists in `reference_sequences`, this returns `(<index in the list>,
    /// <the matched reference sequence>, <start position>, <end position>)`; otherwise, `None`.
    ///
    /// The start and end positions are assumed to be 1-based.
    pub fn resolve<'a>(
        &self,
        reference_sequences: &'a ReferenceSequences,
    ) -> Option<(usize, &'a ReferenceSequence, u64, u64)> {
        match self {
            Self::Mapped { name, start, end } => {
                let (i, _, reference_sequence) = reference_sequences.get_full(name).unwrap();
                let resolved_end = end.unwrap_or(reference_sequence.len() as u64);
                Some((i, reference_sequence, *start, resolved_end))
            }
            Self::Unmapped | Self::All => None,
        }
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Mapped { name, start, end } => {
                write!(f, "{}:{}", name, start)?;

                if let Some(e) = end {
                    write!(f, "-{}", e)?;
                }

                Ok(())
            }
            Self::Unmapped => write!(f, "{}", UNMAPPED_NAME),
            Self::All => write!(f, "{}", ALL_NAME),
        }
    }
}

#[derive(Debug)]
pub enum ParseError {
    MissingReferenceSequenceName,
    InvalidStartPosition(Box<dyn std::error::Error>),
    InvalidEndPosition(Box<dyn std::error::Error>),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => write!(f, "invalid region"),
            Self::InvalidStartPosition(e) => write!(f, "invalid start position: {}", e),
            Self::InvalidEndPosition(e) => write!(f, "invalid end position: {}", e),
        }
    }
}

impl FromStr for Region {
    type Err = ParseError;

    /// Parses a string to a region.
    ///
    /// A region string is specified as
    /// `<reference-sequence-name>[:<start-position>[-<end-position>]]`.
    ///
    /// The reference sequence name can be "*" to represent unmapped records; or ".", all records.
    /// Otherwise, the reference sequence name is assumed to exist in the reference sequence
    /// dictionary.
    ///
    /// If no start position is given, the minimum position of 1 is used. If no end position is
    /// given, the region is unbounded, and the maximum position of the reference sequence, i.e.,
    /// its length, is expected to be used.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == UNMAPPED_NAME {
            return Ok(Self::Unmapped);
        } else if s == ALL_NAME {
            return Ok(Self::All);
        }

        let mut components = s.split(|c| c == ':' || c == '-');

        let name = components
            .next()
            .map(|t| t.into())
            .ok_or_else(|| ParseError::MissingReferenceSequenceName)?;

        let start = match components.next() {
            Some(t) => t
                .parse()
                .map_err(|e| ParseError::InvalidStartPosition(Box::new(e)))?,
            None => MIN_POSITION,
        };

        let end = match components.next() {
            Some(t) => t
                .parse()
                .map(Some)
                .map_err(|e| ParseError::InvalidEndPosition(Box::new(e)))?,
            None => None,
        };

        Ok(Self::Mapped { name, start, end })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve() {
        let reference_sequences: ReferenceSequences = vec![
            (
                String::from("sq0"),
                ReferenceSequence::new(String::from("sq0"), 8),
            ),
            (
                String::from("sq1"),
                ReferenceSequence::new(String::from("sq1"), 13),
            ),
            (
                String::from("sq0"),
                ReferenceSequence::new(String::from("sq2"), 21),
            ),
        ]
        .into_iter()
        .collect();

        let region = Region::mapped("sq1", 5, Some(8));
        let actual = region.resolve(&reference_sequences);
        let expected = Some((1, &reference_sequences["sq1"], 5, 8));
        assert_eq!(actual, expected);

        let region = Region::mapped("sq1", 5, None);
        let actual = region.resolve(&reference_sequences);
        let expected = Some((1, &reference_sequences["sq1"], 5, 13));
        assert_eq!(actual, expected);

        let region = Region::Unmapped;
        assert_eq!(region.resolve(&reference_sequences), None);

        let region = Region::All;
        assert_eq!(region.resolve(&reference_sequences), None);
    }

    #[test]
    fn test_fmt() {
        let region = Region::mapped("sq2", 3, Some(5));
        assert_eq!(format!("{}", region), "sq2:3-5");

        let region = Region::mapped("sq2", 3, None);
        assert_eq!(format!("{}", region), "sq2:3");

        assert_eq!(format!("{}", Region::Unmapped), "*");

        assert_eq!(format!("{}", Region::All), ".");
    }

    #[test]
    fn test_from_str_reference_sequences() {
        match "sq2:3-5".parse() {
            Ok(Region::Mapped { name, start, end }) => {
                assert_eq!(name, "sq2");
                assert_eq!(start, 3);
                assert_eq!(end, Some(5));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_from_str_reference_sequences_with_no_end() {
        match "sq2:3".parse() {
            Ok(Region::Mapped { name, start, end }) => {
                assert_eq!(name, "sq2");
                assert_eq!(start, 3);
                assert_eq!(end, None);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_from_str_reference_sequences_with_no_start_or_end() {
        match "sq2".parse() {
            Ok(Region::Mapped { name, start, end }) => {
                assert_eq!(name, "sq2");
                assert_eq!(start, 1);
                assert_eq!(end, None);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_from_str_for_unmapped_region() {
        let region = "*".parse().unwrap();
        assert!(matches!(region, Region::Unmapped));
    }

    #[test]
    fn test_from_str_for_all_regions() {
        let region = ".".parse().unwrap();
        assert!(matches!(region, Region::All));
    }
}
