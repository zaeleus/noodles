use std::io;

use super::header::{ReferenceSequence, ReferenceSequences};

/// SAM(-like) record extensions.
pub trait RecordExt {
    /// Returns the associated reference sequence.
    fn reference_sequence<'rs>(
        &self,
        reference_sequences: &'rs ReferenceSequences,
    ) -> Option<io::Result<&'rs ReferenceSequence>>;
}
