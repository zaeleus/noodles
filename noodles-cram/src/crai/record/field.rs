#[derive(Clone, Copy, Debug)]
pub enum Field {
    ReferenceSequenceId,
    AlignmentStart,
    AlignmentSpan,
    Offset,
    Landmark,
    SliceLength,
}
