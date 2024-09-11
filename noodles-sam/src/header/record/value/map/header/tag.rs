//! SAM header header tag.

use std::marker::PhantomData;

use crate::header::record::value::map::{
    self,
    tag::{Other, LENGTH},
};

/// A header tag.
pub type Tag = map::tag::Tag<Standard>;

pub(crate) const VERSION: Tag = map::tag::Tag::Standard(Standard::Version);

/// Sort order (`SO`).
pub const SORT_ORDER: Other<Standard> = Other([b'S', b'O'], PhantomData);

/// Group order (`GO`).
pub const GROUP_ORDER: Other<Standard> = Other([b'G', b'O'], PhantomData);

/// Subsort order (`SS`).
pub const SUBSORT_ORDER: Other<Standard> = Other([b'S', b'S'], PhantomData);

const VN: [u8; LENGTH] = [b'V', b'N'];

/// A SAM header header tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Format version (`VN`).
    Version,
}

impl map::tag::Standard for Standard {}

impl AsRef<[u8; LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::Version => &VN,
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match b {
            VN => Ok(Self::Version),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Version => VN,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Version.as_ref(), b"VN");
    }

    #[test]
    fn test_try_from_u8_array_for_standard() {
        assert_eq!(Standard::try_from([b'V', b'N']), Ok(Standard::Version));
        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; LENGTH]>::from(Standard::Version), [b'V', b'N']);
    }
}
