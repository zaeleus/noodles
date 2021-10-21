//! BAM record reference sequence ID.

use std::{error, fmt};

/// The raw unmapped reference sequence ID.
pub const UNMAPPED: i32 = -1;

const MIN: i32 = 0;

/// A BAM record reference sequence ID.
///
/// A reference sequence ID is the the index of the associated reference sequence in the reference
/// sequence dictionary.
///
/// A value of -1 is used for an unmapped record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceSequenceId(i32);

/// An error returned when a raw SAM record position fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromIntError(i32);

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid value: {}", self.0)
    }
}

impl TryFrom<i32> for ReferenceSequenceId {
    type Error = TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        if n < MIN {
            Err(TryFromIntError(n))
        } else {
            Ok(Self(n))
        }
    }
}

impl From<ReferenceSequenceId> for i32 {
    fn from(reference_sequence_id: ReferenceSequenceId) -> Self {
        reference_sequence_id.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_i32_for_reference_sequence_id() {
        assert_eq!(ReferenceSequenceId::try_from(0), Ok(ReferenceSequenceId(0)));
        assert_eq!(
            ReferenceSequenceId::try_from(13),
            Ok(ReferenceSequenceId(13))
        );

        assert_eq!(ReferenceSequenceId::try_from(-1), Err(TryFromIntError(-1)));
    }

    #[test]
    fn test_from_reference_sequence_id_for_i32() {
        assert_eq!(i32::from(ReferenceSequenceId(0)), 0);
        assert_eq!(i32::from(ReferenceSequenceId(13)), 13);
    }
}
