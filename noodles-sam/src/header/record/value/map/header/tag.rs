//! SAM header header tag.

use crate::header::record::value::map::tag::{self, LENGTH};

/// A SAM header header tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Format version (`VN`).
    Version,
    /// Sorting order of alignments (`SO`).
    SortOrder,
    /// Group order of alignments (`GO`).
    GroupOrder,
    /// Subsort order of alignments (`SS`).
    SubsortOrder,
}

impl tag::Standard for Standard {}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"VN" => Ok(Self::Version),
            b"SO" => Ok(Self::SortOrder),
            b"GO" => Ok(Self::GroupOrder),
            b"SS" => Ok(Self::SubsortOrder),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_array_for_standard() {
        assert_eq!(Standard::try_from([b'V', b'N']), Ok(Standard::Version));
        assert_eq!(Standard::try_from([b'S', b'O']), Ok(Standard::SortOrder));
        assert_eq!(Standard::try_from([b'G', b'O']), Ok(Standard::GroupOrder));
        assert_eq!(Standard::try_from([b'S', b'S']), Ok(Standard::SubsortOrder));

        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }
}
