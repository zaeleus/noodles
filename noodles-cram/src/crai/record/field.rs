#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    ReferenceSequenceId,
    AlignmentStart,
    AlignmentSpan,
    Offset,
    Landmark,
    SliceLength,
}
