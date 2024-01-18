use std::io;

use noodles_sam::{self as sam, alignment::record::fields::ReferenceSequenceId as _};

/// A raw BAM record reference sequence ID.
#[derive(Debug, Eq, PartialEq)]
pub struct ReferenceSequenceId(i32);

impl ReferenceSequenceId {
    pub(super) fn new(n: i32) -> Self {
        Self(n)
    }
}

impl sam::alignment::record::fields::ReferenceSequenceId for ReferenceSequenceId {
    fn try_to_usize(&self) -> io::Result<usize> {
        usize::try_from(self.0).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
        reference_sequence_id.try_to_usize()
    }
}
