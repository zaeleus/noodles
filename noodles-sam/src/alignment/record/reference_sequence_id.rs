use std::io;

/// An alignment record reference sequence ID.
pub trait ReferenceSequenceId {
    /// Converts a raw reference sequence ID to a `usize`.
    fn try_to_usize(&self) -> io::Result<usize>;
}

impl TryFrom<&dyn ReferenceSequenceId> for usize {
    type Error = io::Error;

    fn try_from(raw_reference_sequence_id: &dyn ReferenceSequenceId) -> Result<Self, Self::Error> {
        raw_reference_sequence_id.try_to_usize()
    }
}
