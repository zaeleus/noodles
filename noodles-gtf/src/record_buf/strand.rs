/// A GTF record strand.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    /// Forward strand (`+`).
    Forward,
    /// Reverse strand (`-`).
    Reverse,
}
