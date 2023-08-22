use std::str::FromStr;

use super::ParseError;

/// A VCF header info key.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Standard {
    /// Ancestral allele (`AA`).
    AncestralAllele,
    /// Allele count in genotypes, for each ALT allele, in the same order as listed (`AC`).
    AlleleCount,
    /// Total read depth for each allele (`AD`).
    TotalReadDepths,
    /// Read depth for each allele on the forward strand (`ADF`).
    ForwardStrandReadDepths,
    /// Read depth for each allele on the reverse strand (`ADR`).
    ReverseStrandReadDepths,
    /// Allele frequency for each ALT allele in the same order as listed (`AF`).
    AlleleFrequencies,
    /// Total number of alleles in called genotypes (`AN`).
    TotalAlleleCount,
    /// RMS base quality (`BQ`).
    BaseQuality,
    /// Cigar string describing how to align an alternate allele to the reference allele (`CIGAR`).
    Cigar,
    /// dbSNP membership (`DB`).
    IsInDbSnp,
    /// Combined depth across samples (`DP`).
    TotalDepth,
    // /// End position on CHROM (`END`).
    // EndPosition,
    /// HapMap2 membership (`H2`).
    IsInHapMap2,
    /// HapMap3 membership (`H3`).
    IsInHapMap3,
    /// RMS mapping quality (`MQ`).
    MappingQuality,
    /// Number of MAPQ == 0 reads (`MQ0`).
    ZeroMappingQualityCount,
    /// Number of samples with data (`NS`).
    SamplesWithDataCount,
    /// Strand bias (`SB`).
    StrandBias,
    /// Somatic mutation (`SOMATIC`).
    IsSomaticMutation,
    /// Validated by follow-up experiment (`VALIDATED`).
    IsValidated,
    /// 1000 Genomes membership (`1000G`).
    IsIn1000Genomes,

    /// Imprecise structural variation (`IMPRECISE`).
    IsImprecise,
    /// Indicates a novel structural variation (`NOVEL`).
    IsNovel,
    /// End position of the variant described in this record (`END`).
    EndPosition,
    /// Type of structural variant (`SVTYPE`).
    ///
    /// Deprecated in VCF 4.4.
    SvType,
    /// Difference in length between REF and ALT alleles (`SVLEN`).
    SvLengths,
    /// Confidence interval around POS for imprecise variants (`CIPOS`).
    PositionConfidenceIntervals,
    /// Confidence interval around END for imprecise variants (`CIEND`).
    EndConfidenceIntervals,
    /// Length of base pair identical micro-homology at event breakpoints (`HOMLEN`).
    MicrohomologyLengths,
    /// Sequence of base pair identical micro-homology at event breakpoints (`HOMSEQ`).
    MicrohomologySequences,
    /// ID of the assembled alternate allele in the assembly file (`BKPTID`).
    BreakpointIds,
    /// Mobile element info of the form NAME,START,END,POLARITY (`MEINFO`).
    MobileElementInfo,
    /// Mobile element transduction info of the form CHR,START,END,POLARITY (`METRANS`).
    MobileElementTransductionInfo,
    /// ID of this element in Database of Genomic Variation (`DBVID`).
    DbvId,
    /// ID of this element in DBVAR (`DBVARID`).
    DbVarId,
    /// ID of this element in DBRIP (`DBRIPID`).
    DbRipId,
    /// ID of mate breakends (`MATEID`).
    MateBreakendIds,
    /// ID of partner breakend (`PARID`).
    PartnerBreakendId,
    /// ID of event associated to breakend (`EVENT`).
    BreakendEventId,
    /// ID of associated event (`EVENTTYPE`).
    ///
    /// Added in VCF 4.4.
    EventType,
    /// Confidence interval around the inserted material between breakends (`CILEN`).
    BreakendConfidenceIntervals,
    // /// Read Depth of segment containing breakend (`DP`).
    // BreakendReadDepth,
    /// Read Depth of adjacency (`DPADJ`).
    ///
    /// Removed in VCF 4.4.
    AdjacentReadDepths,
    /// Copy number of segment containing breakend (`CN`).
    BreakendCopyNumber,
    /// Copy number of adjacency (`CNADJ`).
    ///
    /// Removed in VCF 4.4.
    AdjacentCopyNumber,
    /// Confidence interval around copy number for the segment (`CICN`).
    CopyNumberConfidenceIntervals,
    /// Confidence interval around copy number for the adjacency (`CICNADJ`).
    ///
    /// Removed in VCF 4.4.
    AdjacentCopyNumberConfidenceIntervals,
    /// Claim made by the structural variant call. Valid values are D, J, DJ for abundance,
    /// adjacency and both respectively.
    ///
    /// Added in VCF 4.4.
    SvClaim,
    /// Total number of repeat sequences in this allele (`RN`).
    ///
    /// Added in VCF 4.4.
    TotalRepeatSequenceCounts,
    /// Repeat unit sequence of the corresponding repeat sequence (`RUS`).
    ///
    /// Added in VCF 4.4.
    RepeatUnitSequences,
    /// Repeat unit length of the corresponding repeat sequence (`RUL`).
    ///
    /// Added in VCF 4.4.
    RepeatUnitLengths,
    /// Repeat unit count of corresponding repeat sequence (`RUC`).
    ///
    /// Added in VCF 4.4.
    RepeatUnitCounts,
    /// Total number of bases in the corresponding repeat sequence (`RB`).
    ///
    /// Added in VCF 4.4.
    TotalRepeatSequenceBaseCounts,
    /// Confidence interval around RUC (`CIRUC`).
    ///
    /// Added in VCF 4.4.
    RepeatUnitCountConfidenceIntervals,
    /// Confidence interval around RB (`CIRB`).
    ///
    /// Added in VCF 4.4.
    TotalRepeatSequenceBaseCountConfidenceIntervals,
    /// Number of bases in each individual repeat unit (`RUB`).
    ///
    /// Added in VCF 4.4.
    RepeatUnitBaseCounts,
}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::AncestralAllele => "AA",
            Self::AlleleCount => "AC",
            Self::TotalReadDepths => "AD",
            Self::ForwardStrandReadDepths => "ADF",
            Self::ReverseStrandReadDepths => "ADR",
            Self::AlleleFrequencies => "AF",
            Self::TotalAlleleCount => "AN",
            Self::BaseQuality => "BQ",
            Self::Cigar => "CIGAR",
            Self::IsInDbSnp => "DB",
            Self::TotalDepth => "DP",
            // Self::EndPosition => "END",
            Self::IsInHapMap2 => "H2",
            Self::IsInHapMap3 => "H3",
            Self::MappingQuality => "MQ",
            Self::ZeroMappingQualityCount => "MQ0",
            Self::SamplesWithDataCount => "NS",
            Self::StrandBias => "SB",
            Self::IsSomaticMutation => "SOMATIC",
            Self::IsValidated => "VALIDATED",
            Self::IsIn1000Genomes => "1000G",

            Self::IsImprecise => "IMPRECISE",
            Self::IsNovel => "NOVEL",
            Self::EndPosition => "END",
            Self::SvType => "SVTYPE",
            Self::SvLengths => "SVLEN",
            Self::PositionConfidenceIntervals => "CIPOS",
            Self::EndConfidenceIntervals => "CIEND",
            Self::MicrohomologyLengths => "HOMLEN",
            Self::MicrohomologySequences => "HOMSEQ",
            Self::BreakpointIds => "BKPTID",
            Self::MobileElementInfo => "MEINFO",
            Self::MobileElementTransductionInfo => "METRANS",
            Self::DbvId => "DGVID",
            Self::DbVarId => "DBVARID",
            Self::DbRipId => "DBRIPID",
            Self::MateBreakendIds => "MATEID",
            Self::PartnerBreakendId => "PARID",
            Self::BreakendEventId => "EVENT",
            Self::EventType => "EVENTTYPE",
            Self::BreakendConfidenceIntervals => "CILEN",
            // Self::BreakendReadDepth => "DP",
            Self::AdjacentReadDepths => "DPADJ",
            Self::BreakendCopyNumber => "CN",
            Self::AdjacentCopyNumber => "CNADJ",
            Self::CopyNumberConfidenceIntervals => "CICN",
            Self::AdjacentCopyNumberConfidenceIntervals => "CICNADJ",
            Self::SvClaim => "SVCLAIM",
            Self::TotalRepeatSequenceCounts => "RN",
            Self::RepeatUnitSequences => "RUS",
            Self::RepeatUnitLengths => "RUL",
            Self::RepeatUnitCounts => "RUC",
            Self::TotalRepeatSequenceBaseCounts => "RB",
            Self::RepeatUnitCountConfidenceIntervals => "CIRUC",
            Self::TotalRepeatSequenceBaseCountConfidenceIntervals => "CIRB",
            Self::RepeatUnitBaseCounts => "RUB",
        }
    }
}

impl FromStr for Standard {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "AA" => Ok(Self::AncestralAllele),
            "AC" => Ok(Self::AlleleCount),
            "AD" => Ok(Self::TotalReadDepths),
            "ADF" => Ok(Self::ForwardStrandReadDepths),
            "ADR" => Ok(Self::ReverseStrandReadDepths),
            "AF" => Ok(Self::AlleleFrequencies),
            "AN" => Ok(Self::TotalAlleleCount),
            "BQ" => Ok(Self::BaseQuality),
            "CIGAR" => Ok(Self::Cigar),
            "DB" => Ok(Self::IsInDbSnp),
            "DP" => Ok(Self::TotalDepth),
            // "END" => Ok(Self::EndPosition),
            "H2" => Ok(Self::IsInHapMap2),
            "H3" => Ok(Self::IsInHapMap3),
            "MQ" => Ok(Self::MappingQuality),
            "MQ0" => Ok(Self::ZeroMappingQualityCount),
            "NS" => Ok(Self::SamplesWithDataCount),
            "SB" => Ok(Self::StrandBias),
            "SOMATIC" => Ok(Self::IsSomaticMutation),
            "VALIDATED" => Ok(Self::IsValidated),
            "1000G" => Ok(Self::IsIn1000Genomes),

            "IMPRECISE" => Ok(Self::IsImprecise),
            "NOVEL" => Ok(Self::IsNovel),
            "END" => Ok(Self::EndPosition),
            "SVTYPE" => Ok(Self::SvType),
            "SVLEN" => Ok(Self::SvLengths),
            "CIPOS" => Ok(Self::PositionConfidenceIntervals),
            "CIEND" => Ok(Self::EndConfidenceIntervals),
            "HOMLEN" => Ok(Self::MicrohomologyLengths),
            "HOMSEQ" => Ok(Self::MicrohomologySequences),
            "BKPTID" => Ok(Self::BreakpointIds),
            "MEINFO" => Ok(Self::MobileElementInfo),
            "METRANS" => Ok(Self::MobileElementTransductionInfo),
            "DGVID" => Ok(Self::DbvId),
            "DBVARID" => Ok(Self::DbVarId),
            "DBRIPID" => Ok(Self::DbRipId),
            "MATEID" => Ok(Self::MateBreakendIds),
            "PARID" => Ok(Self::PartnerBreakendId),
            "EVENT" => Ok(Self::BreakendEventId),
            "EVENTTYPE" => Ok(Self::EventType),
            "CILEN" => Ok(Self::BreakendConfidenceIntervals),
            // "DP" => Ok(Self::BreakendReadDepth),
            "DPADJ" => Ok(Self::AdjacentReadDepths),
            "CN" => Ok(Self::BreakendCopyNumber),
            "CNADJ" => Ok(Self::AdjacentCopyNumber),
            "CICN" => Ok(Self::CopyNumberConfidenceIntervals),
            "CICNADJ" => Ok(Self::AdjacentCopyNumberConfidenceIntervals),
            "SVCLAIM" => Ok(Self::SvClaim),
            "RN" => Ok(Self::TotalRepeatSequenceCounts),
            "RUS" => Ok(Self::RepeatUnitSequences),
            "RUL" => Ok(Self::RepeatUnitLengths),
            "RUC" => Ok(Self::RepeatUnitCounts),
            "RB" => Ok(Self::TotalRepeatSequenceBaseCounts),
            "CIRUC" => Ok(Self::RepeatUnitCountConfidenceIntervals),
            "CIRB" => Ok(Self::TotalRepeatSequenceBaseCountConfidenceIntervals),
            "RUB" => Ok(Self::RepeatUnitBaseCounts),

            _ => Err(ParseError::Invalid),
        }
    }
}
