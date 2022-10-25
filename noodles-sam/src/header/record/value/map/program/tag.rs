//! SAM header program tag.

use crate::header::record::value::map::tag::{self, LENGTH};

const ID: [u8; LENGTH] = [b'I', b'D'];
const PN: [u8; LENGTH] = [b'P', b'N'];
const CL: [u8; LENGTH] = [b'C', b'L'];
const PP: [u8; LENGTH] = [b'P', b'P'];
const DS: [u8; LENGTH] = [b'D', b'S'];
const VN: [u8; LENGTH] = [b'V', b'N'];

/// A SAM header program tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Program ID (`ID`).
    Id,
    /// Program name (`PN`).
    Name,
    /// Command line (`CL`).
    CommandLine,
    /// Previous program ID (`PP`).
    PreviousId,
    /// Description (`DS`).
    Description,
    /// Program version (`VN`).
    Version,
}

impl tag::Standard for Standard {}

impl AsRef<[u8; LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Standard::Id => &ID,
            Standard::Name => &PN,
            Standard::CommandLine => &CL,
            Standard::PreviousId => &PP,
            Standard::Description => &DS,
            Standard::Version => &VN,
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match b {
            ID => Ok(Self::Id),
            PN => Ok(Self::Name),
            CL => Ok(Self::CommandLine),
            PP => Ok(Self::PreviousId),
            DS => Ok(Self::Description),
            VN => Ok(Self::Version),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Id => ID,
            Standard::Name => PN,
            Standard::CommandLine => CL,
            Standard::PreviousId => PP,
            Standard::Description => DS,
            Standard::Version => VN,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Id.as_ref(), &[b'I', b'D']);
        assert_eq!(Standard::Name.as_ref(), &[b'P', b'N']);
        assert_eq!(Standard::CommandLine.as_ref(), &[b'C', b'L']);
        assert_eq!(Standard::PreviousId.as_ref(), &[b'P', b'P']);
        assert_eq!(Standard::Description.as_ref(), &[b'D', b'S']);
        assert_eq!(Standard::Version.as_ref(), &[b'V', b'N']);
    }

    #[test]
    fn test_try_from_u8_array_for_standard() {
        assert_eq!(Standard::try_from([b'I', b'D']), Ok(Standard::Id));
        assert_eq!(Standard::try_from([b'P', b'N']), Ok(Standard::Name));
        assert_eq!(Standard::try_from([b'C', b'L']), Ok(Standard::CommandLine));
        assert_eq!(Standard::try_from([b'P', b'P']), Ok(Standard::PreviousId));
        assert_eq!(Standard::try_from([b'D', b'S']), Ok(Standard::Description));
        assert_eq!(Standard::try_from([b'V', b'N']), Ok(Standard::Version));

        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; LENGTH]>::from(Standard::Id), [b'I', b'D']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Name), [b'P', b'N']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::CommandLine), [b'C', b'L']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::PreviousId), [b'P', b'P']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Description), [b'D', b'S']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Version), [b'V', b'N']);
    }
}
