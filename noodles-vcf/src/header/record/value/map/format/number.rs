/// A VCF header format record number value.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Number {
    /// An explicit size.
    Count(usize),
    /// The number of alternate bases (`A`).
    AlternateBases,
    /// The number of reference and alternate bases (`R`).
    ReferenceAlternateBases,
    /// The number of samples (`G`).
    Samples,
    /// The size is unknown.
    Unknown,
}

impl Default for Number {
    fn default() -> Self {
        Self::Count(1)
    }
}
