//! SAM header reference sequence tag.

use crate::header::record::value::map::tag::{self, LENGTH};

/// A SAM header reference sequence tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Standard {
    /// Reference sequence name (`SN`).
    Name,
    /// Reference sequence length (`LN`).
    Length,
    /// Alternate locus (`AH`).
    AlternativeLocus,
    /// Alternate reference sequence names (`AN`).
    AlternativeNames,
    /// Genome assembly ID (`AS`).
    AssemblyId,
    /// Description (`DS`).
    Description,
    /// MD5 checksum of the reference sequence (`M5`).
    Md5Checksum,
    /// Species (`SP`).
    Species,
    /// Molecule topology (`TP`).
    MoleculeTopology,
    /// URI of the reference sequence (`UR`).
    Uri,
}

impl tag::Standard for Standard {}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"SN" => Ok(Self::Name),
            b"LN" => Ok(Self::Length),
            b"AH" => Ok(Self::AlternativeLocus),
            b"AN" => Ok(Self::AlternativeNames),
            b"AS" => Ok(Self::AssemblyId),
            b"DS" => Ok(Self::Description),
            b"M5" => Ok(Self::Md5Checksum),
            b"SP" => Ok(Self::Species),
            b"TP" => Ok(Self::MoleculeTopology),
            b"UR" => Ok(Self::Uri),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(Standard::try_from([b'S', b'N']), Ok(Standard::Name));
        assert_eq!(Standard::try_from([b'L', b'N']), Ok(Standard::Length));
        assert_eq!(
            Standard::try_from([b'A', b'H']),
            Ok(Standard::AlternativeLocus)
        );
        assert_eq!(
            Standard::try_from([b'A', b'N']),
            Ok(Standard::AlternativeNames)
        );
        assert_eq!(Standard::try_from([b'A', b'S']), Ok(Standard::AssemblyId));
        assert_eq!(Standard::try_from([b'D', b'S']), Ok(Standard::Description));
        assert_eq!(Standard::try_from([b'M', b'5']), Ok(Standard::Md5Checksum));
        assert_eq!(Standard::try_from([b'S', b'P']), Ok(Standard::Species));
        assert_eq!(
            Standard::try_from([b'T', b'P']),
            Ok(Standard::MoleculeTopology)
        );
        assert_eq!(Standard::try_from([b'U', b'R']), Ok(Standard::Uri));

        assert_eq!(Standard::try_from([b'N', b'D']), Err(()));
    }
}
