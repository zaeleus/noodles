use std::{error, fmt};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Key {
    RecordsHaveNames,
    AlignmentStartsAreDeltas,
    ExternalReferenceSequenceIsRequired,
    SubstitutionMatrix,
    TagSets,
    // CRAM 4.0: quality score orientation (0 = original/sequencing, 1 = alignment)
    QualityScoreOrientation,
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
            [b'R', b'N'] => Ok(Self::RecordsHaveNames),
            [b'A', b'P'] => Ok(Self::AlignmentStartsAreDeltas),
            [b'R', b'R'] => Ok(Self::ExternalReferenceSequenceIsRequired),
            [b'S', b'M'] => Ok(Self::SubstitutionMatrix),
            [b'T', b'D'] => Ok(Self::TagSets),
            [b'Q', b'O'] => Ok(Self::QualityScoreOrientation),
            _ => Err(TryFromByteArrayError(b)),
        }
    }
}

impl From<Key> for [u8; 2] {
    fn from(key: Key) -> Self {
        match key {
            Key::RecordsHaveNames => [b'R', b'N'],
            Key::AlignmentStartsAreDeltas => [b'A', b'P'],
            Key::ExternalReferenceSequenceIsRequired => [b'R', b'R'],
            Key::SubstitutionMatrix => [b'S', b'M'],
            Key::TagSets => [b'T', b'D'],
            Key::QualityScoreOrientation => [b'Q', b'O'],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_array_for_key() {
        assert_eq!(Key::try_from([b'R', b'N']), Ok(Key::RecordsHaveNames));
        assert_eq!(
            Key::try_from([b'A', b'P']),
            Ok(Key::AlignmentStartsAreDeltas)
        );
        assert_eq!(
            Key::try_from([b'R', b'R']),
            Ok(Key::ExternalReferenceSequenceIsRequired)
        );
        assert_eq!(Key::try_from([b'S', b'M']), Ok(Key::SubstitutionMatrix));
        assert_eq!(Key::try_from([b'T', b'D']), Ok(Key::TagSets));
        assert_eq!(
            Key::try_from([b'Q', b'O']),
            Ok(Key::QualityScoreOrientation)
        );

        assert_eq!(
            Key::try_from([b'Z', b'Z']),
            Err(TryFromByteArrayError([b'Z', b'Z']))
        );
    }

    #[test]
    fn test_from_key_for_u8_array() {
        assert_eq!(<[u8; 2]>::from(Key::RecordsHaveNames), [b'R', b'N']);
        assert_eq!(<[u8; 2]>::from(Key::AlignmentStartsAreDeltas), [b'A', b'P']);
        assert_eq!(
            <[u8; 2]>::from(Key::ExternalReferenceSequenceIsRequired),
            [b'R', b'R']
        );
        assert_eq!(<[u8; 2]>::from(Key::SubstitutionMatrix), [b'S', b'M']);
        assert_eq!(<[u8; 2]>::from(Key::TagSets), [b'T', b'D']);
        assert_eq!(<[u8; 2]>::from(Key::QualityScoreOrientation), [b'Q', b'O']);
    }
}
