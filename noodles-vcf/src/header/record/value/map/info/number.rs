use std::fmt;

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

impl fmt::Display for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Count(n) => write!(f, "{n}"),
            Self::AlternateBases => f.write_str("A"),
            Self::ReferenceAlternateBases => f.write_str("R"),
            Self::Samples => f.write_str("G"),
            Self::Unknown => f.write_str("."),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Number::default(), Number::Count(1));
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Number::Count(1).to_string(), "1");
        assert_eq!(Number::AlternateBases.to_string(), "A");
        assert_eq!(Number::ReferenceAlternateBases.to_string(), "R");
        assert_eq!(Number::Samples.to_string(), "G");
        assert_eq!(Number::Unknown.to_string(), ".");
    }
}
