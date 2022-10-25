//! SAM header header tag.

use crate::header::record::value::map::tag::{self, LENGTH};

const VN: [u8; LENGTH] = [b'V', b'N'];
const SO: [u8; LENGTH] = [b'S', b'O'];
const GO: [u8; LENGTH] = [b'G', b'O'];
const SS: [u8; LENGTH] = [b'S', b'S'];

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

impl AsRef<[u8; LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::Version => &VN,
            Self::SortOrder => &SO,
            Self::GroupOrder => &GO,
            Self::SubsortOrder => &SS,
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match b {
            VN => Ok(Self::Version),
            SO => Ok(Self::SortOrder),
            GO => Ok(Self::GroupOrder),
            SS => Ok(Self::SubsortOrder),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Version => VN,
            Standard::SortOrder => SO,
            Standard::GroupOrder => GO,
            Standard::SubsortOrder => SS,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Version.as_ref(), &[b'V', b'N']);
        assert_eq!(Standard::SortOrder.as_ref(), &[b'S', b'O']);
        assert_eq!(Standard::GroupOrder.as_ref(), &[b'G', b'O']);
        assert_eq!(Standard::SubsortOrder.as_ref(), &[b'S', b'S']);
    }

    #[test]
    fn test_try_from_u8_array_for_standard() {
        assert_eq!(Standard::try_from([b'V', b'N']), Ok(Standard::Version));
        assert_eq!(Standard::try_from([b'S', b'O']), Ok(Standard::SortOrder));
        assert_eq!(Standard::try_from([b'G', b'O']), Ok(Standard::GroupOrder));
        assert_eq!(Standard::try_from([b'S', b'S']), Ok(Standard::SubsortOrder));

        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; LENGTH]>::from(Standard::Version), [b'V', b'N']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::SortOrder), [b'S', b'O']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::GroupOrder), [b'G', b'O']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::SubsortOrder), [b'S', b'S']);
    }
}
