/// A GFF record field.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    /// The reference sequence name.
    ReferenceSequenceName,
    /// The source.
    Source,
    /// The feature type.
    Type,
    /// The start position.
    Start,
    /// The end position.
    End,
    /// The score.
    Score,
    /// The strand.
    Strand,
    /// The phase.
    Phase,
    /// The attributes.
    Attributes,
}
