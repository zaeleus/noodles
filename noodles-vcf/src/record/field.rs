use std::fmt;

#[derive(Clone, Copy, Debug)]
pub enum Field {
    Chromosome,
    Position,
    Id,
    ReferenceBases,
    AlternateBases,
    QualityScore,
    FilterStatus,
    Information,
    Format,
}

impl AsRef<str> for Field {
    fn as_ref(&self) -> &str {
        match self {
            Self::Chromosome => "CHROM",
            Self::Position => "POS",
            Self::Id => "ID",
            Self::ReferenceBases => "REF",
            Self::AlternateBases => "ALT",
            Self::QualityScore => "QUAL",
            Self::FilterStatus => "FILTER",
            Self::Information => "INFO",
            Self::Format => "FORMAT",
        }
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Field::Chromosome.to_string(), "CHROM");
        assert_eq!(Field::Position.to_string(), "POS");
        assert_eq!(Field::Id.to_string(), "ID");
        assert_eq!(Field::ReferenceBases.to_string(), "REF");
        assert_eq!(Field::AlternateBases.to_string(), "ALT");
        assert_eq!(Field::QualityScore.to_string(), "QUAL");
        assert_eq!(Field::FilterStatus.to_string(), "FILTER");
        assert_eq!(Field::Information.to_string(), "INFO");
        assert_eq!(Field::Format.to_string(), "FORMAT");
    }
}
