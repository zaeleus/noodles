/// A VCF number describing the cardinality of a field.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
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

impl Number {
    /// The number of alternate bases (`A`).
    pub const A: Self = Self::AlternateBases;

    /// The number of reference and alternate bases (`R`).
    pub const R: Self = Self::ReferenceAlternateBases;

    /// The number of samples (`G`).
    pub const G: Self = Self::Samples;
}

impl Default for Number {
    fn default() -> Self {
        Self::Count(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Number::default(), Number::Count(1));
    }
}
