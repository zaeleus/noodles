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
    /// The number of local alternate bases (`LA`).
    LocalAlternateBases,
    /// The number of local reference and alternate bases (`LR`).
    LocalReferenceAlternateBases,
    /// The number of local samples (`LG`).
    LocalSamples,
    /// The number of genotype alleles (`P`).
    Ploidy,
    /// The number of base modifications (`M`).
    BaseModifications,
    /// The size is unknown.
    Unknown,
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
