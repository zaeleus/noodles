pub mod reference_sequence;

pub use self::reference_sequence::ReferenceSequence;

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
    /// Creates a coordinate-sorted index (CSI).
    pub fn new(
        min_shift: i32,
        depth: i32,
        aux: Vec<u8>,
        reference_sequences: Vec<ReferenceSequence>,
        n_no_coor: Option<u64>,
    ) -> Self {
        Self {
            min_shift,
            depth,
            aux,
            reference_sequences,
            n_no_coor,
        }
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
        Self {
            min_shift: 14,
            depth: 5,
            aux: Vec::new(),
            reference_sequences: Vec::new(),
            n_no_coor: None,
        }
    }
}
