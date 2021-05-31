//! Coordinate-sorted index and fields.

mod builder;
pub mod reference_sequence;

pub use self::{builder::Builder, reference_sequence::ReferenceSequence};

/// A coordinate-sorted index (CSI).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Index {
    min_shift: i32,
    depth: i32,
    aux: Vec<u8>,
    reference_sequences: Vec<ReferenceSequence>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Returns a builder to create an index from each of its fields.
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the number of bits for the minimum interval.
    pub fn min_shift(&self) -> i32 {
        self.min_shift
    }

    /// Returns the depth of the binning index.
    pub fn depth(&self) -> i32 {
        self.depth
    }

    /// Returns the auxiliary data.
    pub fn aux(&self) -> &[u8] {
        &self.aux
    }

    /// Returns the list of reference sequences.
    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unmapped records in the associated file.
    pub fn unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}

impl Default for Index {
    fn default() -> Self {
        Self::builder().build()
    }
}
