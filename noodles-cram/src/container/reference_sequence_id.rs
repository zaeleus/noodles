use std::num;

const UNMAPPED: i32 = -1;
const MULTIPLE_REFERENCE_SEQUENCES: i32 = -2;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceSequenceId {
    /// A reference sequence ID.
    Some(usize),
    /// Unmapped (-1).
    None,
    /// Multiple reference sequence IDs (-2).
    Many,
}

impl ReferenceSequenceId {
    pub fn is_some(self) -> bool {
        matches!(self, Self::Some(_))
    }

    pub fn is_none(self) -> bool {
        matches!(self, Self::None)
    }

    pub fn is_many(self) -> bool {
        matches!(self, Self::Many)
    }
}

impl Default for ReferenceSequenceId {
    fn default() -> Self {
        Self::None
    }
}

impl TryFrom<i32> for ReferenceSequenceId {
    type Error = num::TryFromIntError;

    fn try_from(n: i32) -> Result<Self, Self::Error> {
        if n == UNMAPPED {
            Ok(Self::None)
        } else if n == MULTIPLE_REFERENCE_SEQUENCES {
            Ok(Self::Many)
        } else {
            usize::try_from(n).map(Self::Some)
        }
    }
}

impl TryFrom<ReferenceSequenceId> for i32 {
    type Error = num::TryFromIntError;

    fn try_from(reference_sequence_id: ReferenceSequenceId) -> Result<Self, Self::Error> {
        match reference_sequence_id {
            ReferenceSequenceId::Some(n) => i32::try_from(n),
            ReferenceSequenceId::None => Ok(UNMAPPED),
            ReferenceSequenceId::Many => Ok(MULTIPLE_REFERENCE_SEQUENCES),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(ReferenceSequenceId::default(), ReferenceSequenceId::None);
    }

    #[test]
    fn test_try_from_i32_for_reference_sequence_id() {
        assert_eq!(
            ReferenceSequenceId::try_from(1),
            Ok(ReferenceSequenceId::Some(1))
        );

        assert_eq!(
            ReferenceSequenceId::try_from(0),
            Ok(ReferenceSequenceId::Some(0))
        );

        assert_eq!(
            ReferenceSequenceId::try_from(-1),
            Ok(ReferenceSequenceId::None)
        );

        assert_eq!(
            ReferenceSequenceId::try_from(-2),
            Ok(ReferenceSequenceId::Many)
        );

        assert!(ReferenceSequenceId::try_from(-3).is_err());
    }

    #[test]
    fn test_try_from_reference_sequence_id_for_i32() {
        assert_eq!(i32::try_from(ReferenceSequenceId::Some(0)), Ok(0));
        assert_eq!(i32::try_from(ReferenceSequenceId::None), Ok(-1));
        assert_eq!(i32::try_from(ReferenceSequenceId::Many), Ok(-2));
    }
}
