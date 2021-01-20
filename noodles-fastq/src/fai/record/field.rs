/// A FASTA index record field.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    /// The read name.
    ReadName,
    /// The total length of the sequence.
    Length,
    /// The offset to the sequence from the start.
    SequenceOffset,
    /// The number of bases in the sequence.
    LineBases,
    /// The number of characters in the sequence.
    LineWidth,
    /// The offset to the quality scores from the start.
    QualityScoresOffset,
}
