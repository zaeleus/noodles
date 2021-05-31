//! Coordinate-sorted index (CSI) reference sequence and fields.

mod bin;

pub use self::bin::Bin;

/// A CSI reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
}

impl ReferenceSequence {
    /// Creates a CSI reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new());
    /// ```
    pub fn new(bins: Vec<Bin>) -> Self {
        Self { bins }
    }

    /// Returns the list of bins in the reference sequence.
    ///
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new());
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }
}
