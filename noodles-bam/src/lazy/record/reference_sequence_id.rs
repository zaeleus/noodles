use std::io;

/// A raw BAM record reference sequence ID.
#[derive(Debug, Eq, PartialEq)]
pub struct ReferenceSequenceId(i32);

impl ReferenceSequenceId {
    pub(super) fn new(n: i32) -> Self {
        Self(n)
    }
}

impl From<ReferenceSequenceId> for i32 {
    fn from(reference_sequence_id: ReferenceSequenceId) -> Self {
        reference_sequence_id.0
    }
}

impl TryFrom<ReferenceSequenceId> for usize {
    type Error = io::Error;

    fn try_from(reference_sequence_id: ReferenceSequenceId) -> Result<Self, Self::Error> {
        usize::try_from(i32::from(reference_sequence_id))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
