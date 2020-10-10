use std::fmt;

/// A SAM record field.
///
/// A SAM record has 11 required fields and 1 optional field (`Data`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Field {
    /// Read name (`QNAME`).
    Name,
    /// Flags (`FLAG`).
    Flags,
    /// Reference sequence name (`RNAME`).
    ReferenceSequenceName,
    /// 1-based start position (`POS`).
    Position,
    /// Mapping quality (`MAPQ`).
    MappingQuality,
    /// CIGAR operations (`CIGAR`).
    Cigar,
    /// Mate reference sequence name (`RNEXT`).
    MateReferenceSequenceName,
    /// 1-based mate start position (`PNEXT`).
    MatePosition,
    /// Template length (`TLEN`).
    TemplateLength,
    /// Sequence (`SEQ`).
    Sequence,
    /// Quality scores (`QUAL`).
    QualityScores,
    /// Optional data.
    Data,
}

impl AsRef<str> for Field {
    fn as_ref(&self) -> &str {
        match self {
            Self::Name => "QNAME",
            Self::Flags => "FLAG",
            Self::ReferenceSequenceName => "RNAME",
            Self::Position => "POS",
            Self::MappingQuality => "MAPQ",
            Self::Cigar => "CIGAR",
            Self::MateReferenceSequenceName => "RNEXT",
            Self::MatePosition => "PNEXT",
            Self::TemplateLength => "TLEN",
            Self::Sequence => "SEQ",
            Self::QualityScores => "QUAL",
            Self::Data => "DATA",
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
        assert_eq!(Field::Name.to_string(), "QNAME");
        assert_eq!(Field::Flags.to_string(), "FLAG");
        assert_eq!(Field::ReferenceSequenceName.to_string(), "RNAME");
        assert_eq!(Field::Position.to_string(), "POS");
        assert_eq!(Field::MappingQuality.to_string(), "MAPQ");
        assert_eq!(Field::Cigar.to_string(), "CIGAR");
        assert_eq!(Field::MateReferenceSequenceName.to_string(), "RNEXT");
        assert_eq!(Field::MatePosition.to_string(), "PNEXT");
        assert_eq!(Field::TemplateLength.to_string(), "TLEN");
        assert_eq!(Field::Sequence.to_string(), "SEQ");
        assert_eq!(Field::QualityScores.to_string(), "QUAL");
        assert_eq!(Field::Data.to_string(), "DATA");
    }
}
