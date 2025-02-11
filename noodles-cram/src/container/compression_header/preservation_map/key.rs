//! Container compression header preservation map key.

use std::{error, fmt};

/// A container compression header preservation map key.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Key {
    /// Read names are preserved for all records (`RN`).
    ReadNamesIncluded,
    /// AP data series is delta (`AP`).
    ApDataSeriesDelta,
    /// A reference sequence is required to restore data (`RR`).
    ReferenceRequired,
    /// Substitution matrix (`SM`).
    SubstitutionMatrix,
    /// A list of tag sets (`TD`).
    TagSets,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteArrayError([u8; 2]);

impl error::Error for TryFromByteArrayError {}

impl fmt::Display for TryFromByteArrayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid preservation map key: {:#x?}", self.0)
    }
}

impl TryFrom<[u8; 2]> for Key {
    type Error = TryFromByteArrayError;

    fn try_from(b: [u8; 2]) -> Result<Self, Self::Error> {
        match b {
            [b'R', b'N'] => Ok(Self::ReadNamesIncluded),
            [b'A', b'P'] => Ok(Self::ApDataSeriesDelta),
            [b'R', b'R'] => Ok(Self::ReferenceRequired),
            [b'S', b'M'] => Ok(Self::SubstitutionMatrix),
            [b'T', b'D'] => Ok(Self::TagSets),
            _ => Err(TryFromByteArrayError(b)),
        }
    }
}

impl From<Key> for [u8; 2] {
    fn from(key: Key) -> Self {
        match key {
            Key::ReadNamesIncluded => [b'R', b'N'],
            Key::ApDataSeriesDelta => [b'A', b'P'],
            Key::ReferenceRequired => [b'R', b'R'],
            Key::SubstitutionMatrix => [b'S', b'M'],
            Key::TagSets => [b'T', b'D'],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_byte_slice_for_key() {
        assert_eq!(Key::try_from([b'R', b'N']), Ok(Key::ReadNamesIncluded));
        assert_eq!(Key::try_from([b'A', b'P']), Ok(Key::ApDataSeriesDelta));
        assert_eq!(Key::try_from([b'R', b'R']), Ok(Key::ReferenceRequired));
        assert_eq!(Key::try_from([b'S', b'M']), Ok(Key::SubstitutionMatrix));
        assert_eq!(Key::try_from([b'T', b'D']), Ok(Key::TagSets));

        assert_eq!(
            Key::try_from([b'Z', b'Z']),
            Err(TryFromByteArrayError([b'Z', b'Z']))
        );
    }

    #[test]
    fn test_from_key_for_u8_array() {
        assert_eq!(<[u8; 2]>::from(Key::ReadNamesIncluded), [b'R', b'N']);
        assert_eq!(<[u8; 2]>::from(Key::ApDataSeriesDelta), [b'A', b'P']);
        assert_eq!(<[u8; 2]>::from(Key::ReferenceRequired), [b'R', b'R']);
        assert_eq!(<[u8; 2]>::from(Key::SubstitutionMatrix), [b'S', b'M']);
        assert_eq!(<[u8; 2]>::from(Key::TagSets), [b'T', b'D']);
    }
}
