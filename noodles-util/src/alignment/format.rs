/// An alignment format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// Sequence Alignment/Map (SAM).
    Sam,
    /// Binary Alignment/Map (BAM).
    Bam,
    /// CRAM.
    Cram,
}
