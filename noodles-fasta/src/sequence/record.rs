//! Sequence record.

/// A sequence record.
pub trait Record {
    /// Returns the name.
    fn name(&self) -> &[u8];

    /// Returns the description.
    fn description(&self) -> Option<&[u8]>;

    /// Returns the sequence.
    fn sequence(&self) -> &[u8];

    /// Returns the quality scores.
    fn quality_scores(&self) -> &[u8];
}
