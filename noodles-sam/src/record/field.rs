#[derive(Clone, Copy, Debug)]
pub enum Field {
    Name,
    Flags,
    ReferenceSequenceName,
    Position,
    MappingQuality,
    Cigar,
    MateReferenceSequenceName,
    MatePosition,
    TemplateLength,
    Sequence,
    QualityScores,
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
