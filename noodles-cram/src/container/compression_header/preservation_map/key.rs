//! Container compression header preservation map key.

use std::{convert::TryFrom, error, fmt};

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
    /// A list of lists of tag IDs (`TD`).
    TagIdsDictionary,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteSliceError;

impl error::Error for TryFromByteSliceError {}

impl fmt::Display for TryFromByteSliceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid key")
    }
}

impl TryFrom<&[u8]> for Key {
    type Error = TryFromByteSliceError;

    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        match bytes {
            b"RN" => Ok(Self::ReadNamesIncluded),
            b"AP" => Ok(Self::ApDataSeriesDelta),
            b"RR" => Ok(Self::ReferenceRequired),
            b"SM" => Ok(Self::SubstitutionMatrix),
            b"TD" => Ok(Self::TagIdsDictionary),
            _ => Err(TryFromByteSliceError),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_byte_slice_for_key() {
        assert_eq!(Key::try_from(&b"RN"[..]), Ok(Key::ReadNamesIncluded));
        assert_eq!(Key::try_from(&b"AP"[..]), Ok(Key::ApDataSeriesDelta));
        assert_eq!(Key::try_from(&b"RR"[..]), Ok(Key::ReferenceRequired));
        assert_eq!(Key::try_from(&b"SM"[..]), Ok(Key::SubstitutionMatrix));
        assert_eq!(Key::try_from(&b"TD"[..]), Ok(Key::TagIdsDictionary));

        assert_eq!(Key::try_from(&b"ZZ"[..]), Err(TryFromByteSliceError));
    }
}
