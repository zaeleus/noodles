use std::{error, fmt, ops::Bound, str::FromStr};

use noodles_sam::header::ReferenceSequence;

// Position coordinates are 1-based.
const MIN_POSITION: u64 = 1;

static UNMAPPED_NAME: &str = "*";
static ALL_NAME: &str = ".";

#[derive(Clone, Debug)]
pub enum Region {
    Mapped {
        name: String,
        start: u64,
        end: Bound<u64>,
    },
    Unmapped,
    All,
}

impl Region {
    /// Creates a new mapped region.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::ops::Bound;
    /// use noodles::Region;
    ///
    /// let name = String::from("sn1");
    /// let region = Region::mapped(name, 1, Bound::Included(5));
    ///
    /// assert!(matches!(region, Region::Mapped {
    ///     name,
    ///     start: 1,
    ///     end: Bound::Included(5),
    /// }));
    /// ```
    pub fn mapped(name: String, start: u64, end: Bound<u64>) -> Region {
        Region::Mapped { name, start, end }
    }

    /// Returns the reference name of the region.
    ///
    /// If the region is unmapped, this returns "*". If the region represents
    /// all records, this returns ".".
    ///
    /// # Examples
    ///
    /// ```
    /// use std::ops::Bound;
    /// use noodles::Region;
    ///
    /// let region = Region::mapped("sn1".into(), 1, Bound::Included(5));
    /// assert_eq!(region.name(), "sn1");
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

    pub fn resolve<'a>(
        &self,
        reference_sequences: &'a [ReferenceSequence],
    ) -> Option<(usize, &'a ReferenceSequence, u64, u64)> {
        match self {
            Self::Mapped { name, start, end } => {
                let (i, reference_sequence) = reference_sequences
                    .iter()
                    .enumerate()
                    .find(|(_, s)| s.name() == name)
                    .unwrap();

                let resolved_end = match end {
                    Bound::Included(e) => *e,
                    Bound::Excluded(_) => unimplemented!(),
                    Bound::Unbounded => reference_sequence.len() as u64,
                };

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

                if let Bound::Included(e) = end {
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
                .map(Bound::Included)
                .map_err(|e| ParseError::InvalidEndPosition(Box::new(e)))?,
            None => Bound::Unbounded,
        };

        Ok(Self::Mapped { name, start, end })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let region = Region::mapped(String::from("sn2"), 3, Bound::Included(5));
        assert_eq!(format!("{}", region), "sn2:3-5");

        let region = Region::mapped(String::from("sn2"), 3, Bound::Unbounded);
        assert_eq!(format!("{}", region), "sn2:3");

        assert_eq!(format!("{}", Region::Unmapped), "*");

        assert_eq!(format!("{}", Region::All), ".");
    }

    #[test]
    fn test_from_str_reference_sequences() {
        match "sn2:3-5".parse() {
            Ok(Region::Mapped { name, start, end }) => {
                assert_eq!(name, "sn2");
                assert_eq!(start, 3);
                assert_eq!(end, Bound::Included(5));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_from_str_reference_sequences_with_no_end() {
        match "sn2:3".parse() {
            Ok(Region::Mapped { name, start, end }) => {
                assert_eq!(name, "sn2");
                assert_eq!(start, 3);
                assert_eq!(end, Bound::Unbounded);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_from_str_reference_sequences_with_no_start_or_end() {
        match "sn2".parse() {
            Ok(Region::Mapped { name, start, end }) => {
                assert_eq!(name, "sn2");
                assert_eq!(start, 1);
                assert_eq!(end, Bound::Unbounded);
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
