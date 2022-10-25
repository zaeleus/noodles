//! SAM header reference sequence tag.

use crate::header::record::value::map::tag::{self, LENGTH};

const SN: [u8; LENGTH] = [b'S', b'N'];
const LN: [u8; LENGTH] = [b'L', b'N'];
const AH: [u8; LENGTH] = [b'A', b'H'];
const AN: [u8; LENGTH] = [b'A', b'N'];
const AS: [u8; LENGTH] = [b'A', b'S'];
const DS: [u8; LENGTH] = [b'D', b'S'];
const M5: [u8; LENGTH] = [b'M', b'5'];
const SP: [u8; LENGTH] = [b'S', b'P'];
const TP: [u8; LENGTH] = [b'T', b'P'];
const UR: [u8; LENGTH] = [b'U', b'R'];

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

impl AsRef<[u8; LENGTH]> for Standard {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Standard::Name => &SN,
            Standard::Length => &LN,
            Standard::AlternativeLocus => &AH,
            Standard::AlternativeNames => &AN,
            Standard::AssemblyId => &AS,
            Standard::Description => &DS,
            Standard::Md5Checksum => &M5,
            Standard::Species => &SP,
            Standard::MoleculeTopology => &TP,
            Standard::Uri => &UR,
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match b {
            SN => Ok(Self::Name),
            LN => Ok(Self::Length),
            AH => Ok(Self::AlternativeLocus),
            AN => Ok(Self::AlternativeNames),
            AS => Ok(Self::AssemblyId),
            DS => Ok(Self::Description),
            M5 => Ok(Self::Md5Checksum),
            SP => Ok(Self::Species),
            TP => Ok(Self::MoleculeTopology),
            UR => Ok(Self::Uri),
            _ => Err(()),
        }
    }
}

impl From<Standard> for [u8; LENGTH] {
    fn from(tag: Standard) -> Self {
        match tag {
            Standard::Name => SN,
            Standard::Length => LN,
            Standard::AlternativeLocus => AH,
            Standard::AlternativeNames => AN,
            Standard::AssemblyId => AS,
            Standard::Description => DS,
            Standard::Md5Checksum => M5,
            Standard::Species => SP,
            Standard::MoleculeTopology => TP,
            Standard::Uri => UR,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_u8_2_array_for_standard() {
        assert_eq!(Standard::Name.as_ref(), &[b'S', b'N']);
        assert_eq!(Standard::Length.as_ref(), &[b'L', b'N']);
        assert_eq!(Standard::AlternativeLocus.as_ref(), &[b'A', b'H']);
        assert_eq!(Standard::AlternativeNames.as_ref(), &[b'A', b'N']);
        assert_eq!(Standard::AssemblyId.as_ref(), &[b'A', b'S']);
        assert_eq!(Standard::Description.as_ref(), &[b'D', b'S']);
        assert_eq!(Standard::Md5Checksum.as_ref(), &[b'M', b'5']);
        assert_eq!(Standard::Species.as_ref(), &[b'S', b'P']);
        assert_eq!(Standard::MoleculeTopology.as_ref(), &[b'T', b'P']);
        assert_eq!(Standard::Uri.as_ref(), &[b'U', b'R']);
    }

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

    #[test]
    fn test_from_standard_for_u8_2_array() {
        assert_eq!(<[u8; LENGTH]>::from(Standard::Name), [b'S', b'N']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Length), [b'L', b'N']);
        assert_eq!(
            <[u8; LENGTH]>::from(Standard::AlternativeLocus),
            [b'A', b'H']
        );
        assert_eq!(
            <[u8; LENGTH]>::from(Standard::AlternativeNames),
            [b'A', b'N']
        );
        assert_eq!(<[u8; LENGTH]>::from(Standard::AssemblyId), [b'A', b'S']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Description), [b'D', b'S']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Md5Checksum), [b'M', b'5']);
        assert_eq!(<[u8; LENGTH]>::from(Standard::Species), [b'S', b'P']);
        assert_eq!(
            <[u8; LENGTH]>::from(Standard::MoleculeTopology),
            [b'T', b'P']
        );
        assert_eq!(<[u8; LENGTH]>::from(Standard::Uri), [b'U', b'R']);
    }
}
