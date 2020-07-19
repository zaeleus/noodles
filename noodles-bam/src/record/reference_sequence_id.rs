use std::ops::Deref;

const UNMAPPED: i32 = -1;
const MIN: i32 = 0;

/// BAM record reference sequence ID.
///
/// A reference sequence ID is the the index of the associated reference sequence in the reference
/// sequence dictionary.
///
/// A value of -1 is used for an unmapped record.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ReferenceSequenceId(Option<i32>);

impl Deref for ReferenceSequenceId {
    type Target = Option<i32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<i32> for ReferenceSequenceId {
    fn from(n: i32) -> Self {
        if n < MIN {
            Self(None)
        } else {
            Self(Some(n))
        }
    }
}

impl From<ReferenceSequenceId> for i32 {
    fn from(reference_sequence_id: ReferenceSequenceId) -> Self {
        match *reference_sequence_id {
            Some(id) => id,
            None => UNMAPPED,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i32_for_reference_sequence_id() {
        assert_eq!(*ReferenceSequenceId::from(-1), None);
        assert_eq!(*ReferenceSequenceId::from(13), Some(13));
    }

    #[test]
    fn test_from_reference_sequence_id_for_i32() {
        assert_eq!(i32::from(ReferenceSequenceId::from(-1)), -1);
        assert_eq!(i32::from(ReferenceSequenceId::from(13)), 13);
    }
}
