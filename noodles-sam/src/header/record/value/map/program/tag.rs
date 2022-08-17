//! SAM header program tag.

use crate::header::record::value::map::tag::{self, LENGTH};

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

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"ID" => Ok(Self::Id),
            b"PN" => Ok(Self::Name),
            b"CL" => Ok(Self::CommandLine),
            b"PP" => Ok(Self::PreviousId),
            b"DS" => Ok(Self::Description),
            b"VN" => Ok(Self::Version),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
