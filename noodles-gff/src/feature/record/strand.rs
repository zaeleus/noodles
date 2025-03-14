/// A feature record strand.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub enum Strand {
    /// Unstranded (`.`).
    #[default]
    None,
    /// Forward strand (`+`).
    Forward,
    /// Reverse strand (`-`).
    Reverse,
    /// Strandedness is relevant but unknown (`?`).
    Unknown,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Strand::default(), Strand::None);
    }
}
