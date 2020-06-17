/// A SAM record field.
///
/// A SAM record has 11 required fields and 1 optional field (`Data`).
#[derive(Clone, Copy, Debug)]
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

impl Field {
    pub fn name(&self) -> &str {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_name() {
        assert_eq!(Field::Name.name(), "QNAME");
        assert_eq!(Field::Flags.name(), "FLAG");
        assert_eq!(Field::ReferenceSequenceName.name(), "RNAME");
        assert_eq!(Field::Position.name(), "POS");
        assert_eq!(Field::MappingQuality.name(), "MAPQ");
        assert_eq!(Field::Cigar.name(), "CIGAR");
        assert_eq!(Field::MateReferenceSequenceName.name(), "RNEXT");
        assert_eq!(Field::MatePosition.name(), "PNEXT");
        assert_eq!(Field::TemplateLength.name(), "TLEN");
        assert_eq!(Field::Sequence.name(), "SEQ");
        assert_eq!(Field::QualityScores.name(), "QUAL");
        assert_eq!(Field::Data.name(), "DATA");
    }
}
