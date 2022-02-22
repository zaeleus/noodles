//! BAM record reference sequence ID.

use std::{error, fmt};

/// The raw unmapped reference sequence ID.
pub const UNMAPPED: i32 = -1;

/// A BAM record reference sequence ID.
///
/// A reference sequence ID is the the index of the associated reference sequence in the reference
/// sequence dictionary.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ReferenceSequenceId(usize);

/// An error returned when a raw SAM record position fails to convert.
#[deprecated(
    since = "0.16.0",
    note = "Convert the raw reference sequence ID to a `usize` first, and use `ReferenceSequenceId::from` instead."
)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromIntError(i32);

#[allow(deprecated)]
impl error::Error for TryFromIntError {}

#[allow(deprecated)]
impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid value: {}", self.0)
    }
}

impl From<usize> for ReferenceSequenceId {
    fn from(n: usize) -> Self {
        Self(n)
    }
}

#[allow(deprecated)]
impl TryFrom<i32> for ReferenceSequenceId {
    type Error = TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        usize::try_from(n).map(Self).map_err(|_| TryFromIntError(n))
    }
}

impl From<ReferenceSequenceId> for i32 {
    #[allow(useless_deprecated)]
    #[deprecated(
        since = "0.16.0",
        note = "Use `i32::try_from(usize::from(reference_sequence_id))` instead."
    )]
    fn from(reference_sequence_id: ReferenceSequenceId) -> Self {
        // FIXME
        reference_sequence_id.0 as i32
    }
}

impl From<ReferenceSequenceId> for usize {
    fn from(reference_sequence_id: ReferenceSequenceId) -> Self {
        reference_sequence_id.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_usize_for_reference_sequence_id() {
        assert_eq!(ReferenceSequenceId::from(0), ReferenceSequenceId(0));
        assert_eq!(ReferenceSequenceId::from(13), ReferenceSequenceId(13));
    }

    #[allow(deprecated)]
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
    fn test_from_reference_sequence_id_for_usize() {
        assert_eq!(usize::from(ReferenceSequenceId(0)), 0);
        assert_eq!(usize::from(ReferenceSequenceId(13)), 13);
    }
}
