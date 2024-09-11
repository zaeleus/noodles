//! SAM header program tag.

use std::marker::PhantomData;

use crate::header::record::value::map::{
    self,
    tag::{Other, LENGTH},
};

pub(crate) type Tag = map::tag::Tag<Standard>;

pub(crate) const ID: Tag = map::tag::Tag::Standard(Standard::Id);

/// Program name (`PN`).
pub const NAME: Other<Standard> = Other([b'P', b'N'], PhantomData);

/// Command line (`CL`).
pub const COMMAND_LINE: Other<Standard> = Other([b'C', b'L'], PhantomData);

/// Previous program ID (`PP`).
pub const PREVIOUS_PROGRAM_ID: Other<Standard> = Other([b'P', b'P'], PhantomData);

/// Description (`DS`).
pub const DESCRIPTION: Other<Standard> = Other([b'D', b'S'], PhantomData);

/// Program version (`VN`).
pub const VERSION: Other<Standard> = Other([b'V', b'N'], PhantomData);

const ID_VALUE: [u8; LENGTH] = [b'I', b'D'];

/// A SAM header program tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Program ID (`ID`).
    Id,
}

impl map::tag::Standard for Standard {}

impl AsRef<[u8; LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Standard::Id => &ID_VALUE,
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match b {
            ID_VALUE => Ok(Self::Id),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Id => ID_VALUE,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Id.as_ref(), b"ID");
    }

    #[test]
    fn test_try_from_u8_array_for_standard() {
        assert_eq!(Standard::try_from([b'I', b'D']), Ok(Standard::Id));
        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; LENGTH]>::from(Standard::Id), [b'I', b'D']);
    }
}
