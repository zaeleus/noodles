/// A FASTA index record field.
#[derive(Clone, Copy, Debug)]
pub enum Field {
    /// Reference sequence name.
    Name,
    /// The total length of the sequence.
    Length,
    /// Offset from the start.
    Offset,
    /// The number of bases in a line.
    LineBases,
    /// The number of characters in a line.
    LineWidth,
}
