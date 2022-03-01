use std::{error, fmt};

const MIN: i32 = 0;
const UNMAPPED: i32 = -1;
const MULTIPLE_REFERENCE_SEQUENCES: i32 = -2;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceSequenceId {
    /// A reference sequence ID.
    Some(i32),
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
        if n >= MIN {
            Ok(Self::Some(n))
        } else if n == UNMAPPED {
            Ok(Self::None)
        } else if n == MULTIPLE_REFERENCE_SEQUENCES {
            Ok(Self::Many)
        } else {
            Err(TryFromIntError(n))
        }
    }
}

impl From<ReferenceSequenceId> for i32 {
    fn from(reference_sequence_id: ReferenceSequenceId) -> Self {
        match reference_sequence_id {
            ReferenceSequenceId::Some(n) => n,
            ReferenceSequenceId::None => UNMAPPED,
            ReferenceSequenceId::Many => MULTIPLE_REFERENCE_SEQUENCES,
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
    fn test_try_from_itf8_for_reference_sequence_id() {
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

        assert_eq!(ReferenceSequenceId::try_from(-3), Err(TryFromIntError(-3)));
    }

    #[test]
    fn test_from_reference_sequence_id_for_i32() {
        assert_eq!(i32::from(ReferenceSequenceId::Some(1)), 1);
        assert_eq!(i32::from(ReferenceSequenceId::Some(0)), 0);
        assert_eq!(i32::from(ReferenceSequenceId::None), -1);
        assert_eq!(i32::from(ReferenceSequenceId::Many), -2);
    }
}
