//! CSI index builder.

use super::{Index, ReferenceSequence};

/// A coordinate-sorted index (CSI) builder.
pub struct Builder {
    min_shift: i32,
    depth: i32,
    aux: Vec<u8>,
    reference_sequences: Vec<ReferenceSequence>,
    n_no_coor: Option<u64>,
}

impl Builder {
    /// Sets a min shift.
    pub fn set_min_shift(mut self, min_shift: i32) -> Self {
        self.min_shift = min_shift;
        self
    }

    /// Sets a max depth.
    pub fn set_depth(mut self, depth: i32) -> Self {
        self.depth = depth;
        self
    }

    /// Set auxiliary data.
    pub fn set_aux(mut self, aux: Vec<u8>) -> Self {
        self.aux = aux;
        self
    }

    /// Sets reference sequences.
    pub fn set_reference_sequences(mut self, reference_sequences: Vec<ReferenceSequence>) -> Self {
        self.reference_sequences = reference_sequences;
        self
    }

    /// Sets an unmapped read count.
    pub fn set_n_no_coor(mut self, n_no_coor: u64) -> Self {
        self.n_no_coor = Some(n_no_coor);
        self
    }

    /// Builds a coordinate-sorted index (CSI).
    pub fn build(self) -> Index {
        Index {
            min_shift: self.min_shift,
            depth: self.depth,
            aux: self.aux,
            reference_sequences: self.reference_sequences,
            n_no_coor: self.n_no_coor,
        }
    }
}

impl Default for Builder {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert_eq!(builder.min_shift, 14);
        assert_eq!(builder.depth, 5);
        assert!(builder.aux.is_empty());
        assert!(builder.reference_sequences.is_empty());
        assert!(builder.n_no_coor.is_none());
    }
}
