//! GFF record phase.

/// A GFF record phase.
///
/// The phase is used for CDS (coding sequence) features to describe where the next codon begins
/// relative to the 5' end.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Phase {
    /// The codon begins at the first nucleotide (`0`).
    Zero,
    /// The codon begins at the second nucleotide (`1`).
    One,
    /// The codon begins at the third nucleotide (`2`).
    Two,
}
