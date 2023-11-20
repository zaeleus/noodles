use super::index::reference_sequence::Metadata;

/// A binning index reference sequence.
pub trait ReferenceSequence {
    /// Returns the optional metadata for the reference sequence.
    fn metadata(&self) -> Option<&Metadata>;
}
