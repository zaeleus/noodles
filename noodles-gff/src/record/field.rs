#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    ReferenceSequenceName,
    Source,
    Feature,
    Start,
    End,
    Score,
    Strand,
    Frame,
    Attributes,
}
