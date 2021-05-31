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
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::builder().set_min_shift(13).build();
    /// assert_eq!(index.min_shift(), 13);
    /// ```
    pub fn set_min_shift(mut self, min_shift: i32) -> Self {
        self.min_shift = min_shift;
        self
    }

    /// Sets a max depth.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::builder().set_depth(8).build();
    /// assert_eq!(index.depth(), 8);
    /// ```
    pub fn set_depth(mut self, depth: i32) -> Self {
        self.depth = depth;
        self
    }

    /// Set auxiliary data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::builder().set_aux(b"ndls".to_vec()).build();
    /// assert_eq!(index.aux(), b"ndls");
    /// ```
    pub fn set_aux(mut self, aux: Vec<u8>) -> Self {
        self.aux = aux;
        self
    }

    /// Sets reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::{self as csi, index::ReferenceSequence};
    ///
    /// let reference_sequences = vec![ReferenceSequence::new(Vec::new())];
    /// let index = csi::Index::builder()
    ///     .set_reference_sequences(reference_sequences.clone())
    ///     .build();
    ///
    /// assert_eq!(index.reference_sequences(), &reference_sequences);
    /// ```
    pub fn set_reference_sequences(mut self, reference_sequences: Vec<ReferenceSequence>) -> Self {
        self.reference_sequences = reference_sequences;
        self
    }

    /// Sets an unmapped read count.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::builder().set_n_no_coor(21).build();
    /// assert_eq!(index.unmapped_read_count(), Some(21));
    /// ```
    pub fn set_n_no_coor(mut self, n_no_coor: u64) -> Self {
        self.n_no_coor = Some(n_no_coor);
        self
    }

    /// Builds a coordinate-sorted index (CSI).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::builder().build();
    /// ```
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
