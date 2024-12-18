/// A GFF record strand.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Strand {
    /// Unstranded (`.`).
    None,
    /// Forward strand (`+`).
    Forward,
    /// Reverse strand (`-`).
    Reverse,
    /// Strandedness is relevant but unknown (`?`).
    Unknown,
}

impl Default for Strand {
    fn default() -> Self {
        Self::None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Strand::default(), Strand::None);
    }
}
