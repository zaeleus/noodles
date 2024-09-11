//! SAM header reference sequence tag.

use std::marker::PhantomData;

use crate::header::record::value::map::{
    self,
    tag::{self, Other},
};

pub(crate) type Tag = map::tag::Tag<Standard>;

pub(crate) const NAME: Tag = map::tag::Tag::Standard(Standard::Name);
pub(crate) const LENGTH: Tag = map::tag::Tag::Standard(Standard::Length);

/// Alternate locus (`AH`).
pub const ALTERNATIVE_LOCUS: Other<Standard> = Other([b'A', b'H'], PhantomData);

/// Alternate reference sequence names (`AN`).
pub const ALTERNATIVE_NAMES: Other<Standard> = Other([b'A', b'N'], PhantomData);

/// Genome assembly ID (`AS`).
pub const ASSEMBLY_ID: Other<Standard> = Other([b'A', b'S'], PhantomData);

/// Description (`DS`).
pub const DESCRIPTION: Other<Standard> = Other([b'D', b'S'], PhantomData);

/// MD5 checksum of the reference sequence (`M5`).
pub const MD5_CHECKSUM: Other<Standard> = Other([b'M', b'5'], PhantomData);

/// Species (`SP`).
pub const SPECIES: Other<Standard> = Other([b'S', b'P'], PhantomData);

/// Molecule topology (`TP`).
pub const MOLECULE_TOPOLOGY: Other<Standard> = Other([b'T', b'P'], PhantomData);

/// URI of the reference sequence (`UR`).
pub const URI: Other<Standard> = Other([b'U', b'R'], PhantomData);

const SN: [u8; tag::LENGTH] = [b'S', b'N'];
const LN: [u8; tag::LENGTH] = [b'L', b'N'];

/// A SAM header reference sequence tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Reference sequence name (`SN`).
    Name,
    /// Reference sequence length (`LN`).
    Length,
}

impl map::tag::Standard for Standard {}

impl AsRef<[u8; tag::LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; tag::LENGTH] {
        match self {
            Standard::Name => &SN,
            Standard::Length => &LN,
        }
    }
}

impl TryFrom<[u8; tag::LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; tag::LENGTH]) -> Result<Self, Self::Error> {
        match b {
            SN => Ok(Self::Name),
            LN => Ok(Self::Length),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; tag::LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Name => SN,
            Standard::Length => LN,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Name.as_ref(), b"SN");
        assert_eq!(Standard::Length.as_ref(), b"LN");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(Standard::try_from([b'S', b'N']), Ok(Standard::Name));
        assert_eq!(Standard::try_from([b'L', b'N']), Ok(Standard::Length));
        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; tag::LENGTH]>::from(Standard::Name), [b'S', b'N']);
        assert_eq!(<[u8; tag::LENGTH]>::from(Standard::Length), [b'L', b'N']);
    }
}
