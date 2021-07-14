/// A CRAM index record field.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    /// Reference sequence ID.
    ReferenceSequenceId,
    /// Alignment start.
    AlignmentStart,
    /// Alignment end.
    AlignmentSpan,
    /// Offset of the container from the start of the stream.
    Offset,
    /// Offset of the slice from the start of the container.
    Landmark,
    /// Size of the slice in bytes.
    SliceLength,
}
