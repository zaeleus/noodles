/// A FASTA index record field.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    /// The name.
    Name,
    /// The total length of the sequence.
    Length,
    /// The offset from the start.
    Offset,
    /// The number of bases in a line.
    LineBases,
    /// The number of characters in a line.
    LineWidth,
}
