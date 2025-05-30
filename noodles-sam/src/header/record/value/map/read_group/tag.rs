//! SAM header read group tag.

use std::marker::PhantomData;

use crate::header::record::value::map::{
    self,
    tag::{LENGTH, Other},
};

pub(crate) type Tag = map::tag::Tag<Standard>;

pub(crate) const ID: Tag = map::tag::Tag::Standard(Standard::Id);

/// Barcode sequence (`BC`).
pub const BARCODE: Other<Standard> = Other([b'B', b'C'], PhantomData);

/// Sequencing center (`CN`).
pub const SEQUENCING_CENTER: Other<Standard> = Other([b'C', b'N'], PhantomData);

/// Description (`DS`).
pub const DESCRIPTION: Other<Standard> = Other([b'D', b'S'], PhantomData);

/// Datetime of run (`DT`).
pub const PRODUCED_AT: Other<Standard> = Other([b'D', b'T'], PhantomData);

/// Flow order (`FO`).
pub const FLOW_ORDER: Other<Standard> = Other([b'F', b'O'], PhantomData);

/// Key sequence (`KS`).
pub const KEY_SEQUENCE: Other<Standard> = Other([b'K', b'S'], PhantomData);

/// Library (`LB`).
pub const LIBRARY: Other<Standard> = Other([b'L', b'B'], PhantomData);

/// Programs used (`PG`).
pub const PROGRAM: Other<Standard> = Other([b'P', b'G'], PhantomData);

/// Predicted median insert size (`PI`).
pub const PREDICTED_MEDIAN_INSERT_SIZE: Other<Standard> = Other([b'P', b'I'], PhantomData);

/// Platform used (`PL`).
pub const PLATFORM: Other<Standard> = Other([b'P', b'L'], PhantomData);

/// Platform model (`PM`).
pub const PLATFORM_MODEL: Other<Standard> = Other([b'P', b'M'], PhantomData);

/// Platform unit (`PU`).
pub const PLATFORM_UNIT: Other<Standard> = Other([b'P', b'U'], PhantomData);

/// Sample (`SM`).
pub const SAMPLE: Other<Standard> = Other([b'S', b'M'], PhantomData);

const ID_VALUE: [u8; LENGTH] = [b'I', b'D'];

/// A SAM header read group tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Read group ID (`ID`).
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
