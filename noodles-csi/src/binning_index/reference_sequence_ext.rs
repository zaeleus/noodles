use noodles_bgzf as bgzf;

use crate::index::reference_sequence::Metadata;

/// An extension that adds methods to binning index reference sequence types.
pub trait ReferenceSequenceExt {
    /// Returns the optional metadata for the reference sequence.
    fn metadata(&self) -> Option<&Metadata>;

    /// Returns the start position of the first record in the last linear bin.
    fn first_record_in_last_linear_bin_start_position(&self) -> Option<bgzf::VirtualPosition>;
}
