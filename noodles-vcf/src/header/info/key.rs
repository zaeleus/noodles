//! VCF header info key.

use super::Type;
use crate::header::Number;

use std::{error, fmt, str::FromStr};

/// A VCF record info field key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Key {
    // § 1.6.1 Fixed Fields (2021-01-13)
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

    // § 3 INFO keys used for structural variants (2021-01-13)
    /// Imprecise structural variation (`IMPRECISE`).
    IsImprecise,
    /// Indicates a novel structural variation (`NOVEL`).
    IsNovel,
    /// End position of the variant described in this record (`END`).
    EndPosition,
    /// Type of structural variant (`SVTYPE`).
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
    /// Confidence interval around the inserted material between breakends (`CILEN`).
    BreakendConfidenceIntervals,
    // /// Read Depth of segment containing breakend (`DP`).
    // BreakendReadDepth,
    /// Read Depth of adjacency (`DPADJ`).
    AdjacentReadDepths,
    /// Copy number of segment containing breakend (`CN`).
    BreakendCopyNumber,
    /// Copy number of adjacency (`CNADJ`).
    AdjacentCopyNumber,
    /// Confidence interval around copy number for the segment (`CICN`).
    CopyNumberConfidenceIntervals,
    /// Confidence interval around copy number for the adjacency (`CICNADJ`).
    AdjacentCopyNumberConfidenceIntervals,

    /// Any other non-reserved key.
    Other(String),
}

impl AsRef<str> for Key {
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
            Self::BreakendConfidenceIntervals => "CILEN",
            // Self::BreakendReadDepth => "DP",
            Self::AdjacentReadDepths => "DPADJ",
            Self::BreakendCopyNumber => "CN",
            Self::AdjacentCopyNumber => "CNADJ",
            Self::CopyNumberConfidenceIntervals => "CICN",
            Self::AdjacentCopyNumberConfidenceIntervals => "CICNADJ",

            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF record info field key fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

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
            "CILEN" => Ok(Self::BreakendConfidenceIntervals),
            // "DP" => Ok(Self::BreakendReadDepth),
            "DPADJ" => Ok(Self::AdjacentReadDepths),
            "CN" => Ok(Self::BreakendCopyNumber),
            "CNADJ" => Ok(Self::AdjacentCopyNumber),
            "CICN" => Ok(Self::CopyNumberConfidenceIntervals),
            "CICNADJ" => Ok(Self::AdjacentCopyNumberConfidenceIntervals),

            _ => {
                if is_valid_name(s) {
                    Ok(Self::Other(s.into()))
                } else {
                    Err(ParseError::Invalid)
                }
            }
        }
    }
}

// § 1.6.1 Fixed fields
fn is_valid_name_char(c: char) -> bool {
    matches!(c, '0'..='9' | 'A'..='Z' | 'a'..='z' | '_' | '.')
}

fn is_valid_name(s: &str) -> bool {
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if !matches!(c, 'A'..='Z' | 'a'..='z' | '_') {
            return false;
        }
    }

    chars.all(is_valid_name_char)
}

pub(super) fn number(key: &Key) -> Option<Number> {
    match key {
        Key::AncestralAllele => Some(Number::Count(1)),
        Key::AlleleCount => Some(Number::A),
        Key::TotalReadDepths => Some(Number::R),
        Key::ForwardStrandReadDepths => Some(Number::R),
        Key::ReverseStrandReadDepths => Some(Number::R),
        Key::AlleleFrequencies => Some(Number::A),
        Key::TotalAlleleCount => Some(Number::Count(1)),
        Key::BaseQuality => Some(Number::Count(1)),
        Key::Cigar => Some(Number::A),
        Key::IsInDbSnp => Some(Number::Count(0)),
        Key::TotalDepth => Some(Number::Count(1)),
        // Key::EndPosition => Some(Number::Count(1)),
        Key::IsInHapMap2 => Some(Number::Count(0)),
        Key::IsInHapMap3 => Some(Number::Count(0)),
        Key::MappingQuality => Some(Number::Count(1)),
        Key::ZeroMappingQualityCount => Some(Number::Count(1)),
        Key::SamplesWithDataCount => Some(Number::Count(1)),
        Key::StrandBias => Some(Number::Count(4)),
        Key::IsSomaticMutation => Some(Number::Count(0)),
        Key::IsValidated => Some(Number::Count(0)),
        Key::IsIn1000Genomes => Some(Number::Count(0)),

        Key::IsImprecise => Some(Number::Count(0)),
        Key::IsNovel => Some(Number::Count(0)),
        Key::EndPosition => Some(Number::Count(1)),
        Key::SvType => Some(Number::Count(1)),
        Key::SvLengths => Some(Number::Unknown),
        Key::PositionConfidenceIntervals => Some(Number::Count(2)),
        Key::EndConfidenceIntervals => Some(Number::Count(2)),
        Key::MicrohomologyLengths => Some(Number::Unknown),
        Key::MicrohomologySequences => Some(Number::Unknown),
        Key::BreakpointIds => Some(Number::Unknown),
        Key::MobileElementInfo => Some(Number::Count(4)),
        Key::MobileElementTransductionInfo => Some(Number::Count(4)),
        Key::DbvId => Some(Number::Count(1)),
        Key::DbVarId => Some(Number::Count(1)),
        Key::DbRipId => Some(Number::Count(1)),
        Key::MateBreakendIds => Some(Number::Unknown),
        Key::PartnerBreakendId => Some(Number::Count(1)),
        Key::BreakendEventId => Some(Number::Count(1)),
        Key::BreakendConfidenceIntervals => Some(Number::Count(2)),
        // Key::BreakendReadDepth => Some(Number::Count(1)),
        Key::AdjacentReadDepths => Some(Number::Unknown),
        Key::BreakendCopyNumber => Some(Number::Count(1)),
        Key::AdjacentCopyNumber => Some(Number::Unknown),
        Key::CopyNumberConfidenceIntervals => Some(Number::Count(2)),
        Key::AdjacentCopyNumberConfidenceIntervals => Some(Number::Unknown),

        Key::Other(_) => None,
    }
}

pub(super) fn ty(key: &Key) -> Option<Type> {
    match key {
        Key::AncestralAllele => Some(Type::String),
        Key::AlleleCount => Some(Type::Integer),
        Key::TotalReadDepths => Some(Type::Integer),
        Key::ForwardStrandReadDepths => Some(Type::Integer),
        Key::ReverseStrandReadDepths => Some(Type::Integer),
        Key::AlleleFrequencies => Some(Type::Float),
        Key::TotalAlleleCount => Some(Type::Integer),
        Key::BaseQuality => Some(Type::Float),
        Key::Cigar => Some(Type::String),
        Key::IsInDbSnp => Some(Type::Flag),
        Key::TotalDepth => Some(Type::Integer),
        // Key::EndPosition => Some(Type::Integer),
        Key::IsInHapMap2 => Some(Type::Flag),
        Key::IsInHapMap3 => Some(Type::Flag),
        Key::MappingQuality => Some(Type::Float),
        Key::ZeroMappingQualityCount => Some(Type::Integer),
        Key::SamplesWithDataCount => Some(Type::Integer),
        Key::StrandBias => Some(Type::Integer),
        Key::IsSomaticMutation => Some(Type::Flag),
        Key::IsValidated => Some(Type::Flag),
        Key::IsIn1000Genomes => Some(Type::Flag),

        Key::IsImprecise => Some(Type::Flag),
        Key::IsNovel => Some(Type::Flag),
        Key::EndPosition => Some(Type::Integer),
        Key::SvType => Some(Type::String),
        Key::SvLengths => Some(Type::Integer),
        Key::PositionConfidenceIntervals => Some(Type::Integer),
        Key::EndConfidenceIntervals => Some(Type::Integer),
        Key::MicrohomologyLengths => Some(Type::Integer),
        Key::MicrohomologySequences => Some(Type::String),
        Key::BreakpointIds => Some(Type::String),
        Key::MobileElementInfo => Some(Type::String),
        Key::MobileElementTransductionInfo => Some(Type::String),
        Key::DbvId => Some(Type::String),
        Key::DbVarId => Some(Type::String),
        Key::DbRipId => Some(Type::String),
        Key::MateBreakendIds => Some(Type::String),
        Key::PartnerBreakendId => Some(Type::String),
        Key::BreakendEventId => Some(Type::String),
        Key::BreakendConfidenceIntervals => Some(Type::Integer),
        // Key::BreakendReadDepth => Some(Type::Integer),
        Key::AdjacentReadDepths => Some(Type::Integer),
        Key::BreakendCopyNumber => Some(Type::Integer),
        Key::AdjacentCopyNumber => Some(Type::Integer),
        Key::CopyNumberConfidenceIntervals => Some(Type::Integer),
        Key::AdjacentCopyNumberConfidenceIntervals => Some(Type::Integer),

        Key::Other(_) => None,
    }
}

pub(super) fn description(key: &Key) -> Option<&str> {
    match key {
        Key::AncestralAllele => Some("Ancestral allele"),
        Key::AlleleCount => {
            Some("Allele count in genotypes, for each ALT allele, in the same order as listed")
        }
        Key::TotalReadDepths => Some("Total read depth for each allele"),
        Key::ForwardStrandReadDepths => Some("Read depth for each allele on the forward strand"),
        Key::ReverseStrandReadDepths => Some("Read depth for each allele on the reverse strand"),
        Key::AlleleFrequencies => {
            Some("Allele frequency for each ALT allele in the same order as listed")
        }
        Key::TotalAlleleCount => Some("Total number of alleles in called genotypes"),
        Key::BaseQuality => Some("RMS base quality"),
        Key::Cigar => {
            Some("Cigar string describing how to align an alternate allele to the reference allele")
        }
        Key::IsInDbSnp => Some("dbSNP membership"),
        Key::TotalDepth => Some("Combined depth across samples"),
        // Key::EndPosition => Some("End position on CHROM"),
        Key::IsInHapMap2 => Some("HapMap2 membership"),
        Key::IsInHapMap3 => Some("HapMap3 membership"),
        Key::MappingQuality => Some("RMS mapping quality"),
        Key::ZeroMappingQualityCount => Some("Number of MAPQ == 0 reads"),
        Key::SamplesWithDataCount => Some("Number of samples with data"),
        Key::StrandBias => Some("Strand bias"),
        Key::IsSomaticMutation => Some("Somatic mutation"),
        Key::IsValidated => Some("Validated by follow-up experiment"),
        Key::IsIn1000Genomes => Some("1000 Genomes membership"),

        Key::IsImprecise => Some("Imprecise structural variation"),
        Key::IsNovel => Some("Indicates a novel structural variation"),
        Key::EndPosition => Some("End position of the variant described in this record"),
        Key::SvType => Some("Type of structural variant"),
        Key::SvLengths => Some("Difference in length between REF and ALT alleles"),
        Key::PositionConfidenceIntervals => {
            Some("Confidence interval around POS for imprecise variants")
        }
        Key::EndConfidenceIntervals => {
            Some("Confidence interval around END for imprecise variants")
        }
        Key::MicrohomologyLengths => {
            Some("Length of base pair identical micro-homology at event breakpoints")
        }
        Key::MicrohomologySequences => {
            Some("Sequence of base pair identical micro-homology at event breakpoints")
        }
        Key::BreakpointIds => Some("ID of the assembled alternate allele in the assembly file"),
        Key::MobileElementInfo => Some("Mobile element info of the form NAME,START,END,POLARITY"),
        Key::MobileElementTransductionInfo => {
            Some("Mobile element transduction info of the form CHR,START,END,POLARITY")
        }
        Key::DbvId => Some("ID of this element in Database of Genomic Variation"),
        Key::DbVarId => Some("ID of this element in DBVAR"),
        Key::DbRipId => Some("ID of this element in DBRIP"),
        Key::MateBreakendIds => Some("ID of mate breakends"),
        Key::PartnerBreakendId => Some("ID of partner breakend"),
        Key::BreakendEventId => Some("ID of event associated to breakend"),
        Key::BreakendConfidenceIntervals => {
            Some("Confidence interval around the inserted material between breakends")
        }
        // Key::BreakendReadDepth => Some("Read Depth of segment containing breakend"),
        Key::AdjacentReadDepths => Some("Read Depth of adjacency"),
        Key::BreakendCopyNumber => Some("Copy number of segment containing breakend"),
        Key::AdjacentCopyNumber => Some("Copy number of adjacency"),
        Key::CopyNumberConfidenceIntervals => {
            Some("Confidence interval around copy number for the segment")
        }
        Key::AdjacentCopyNumberConfidenceIntervals => {
            Some("Confidence interval around copy number for the adjacency")
        }

        Key::Other(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::AncestralAllele.to_string(), "AA");
        assert_eq!(Key::AlleleCount.to_string(), "AC");
        assert_eq!(Key::TotalReadDepths.to_string(), "AD");
        assert_eq!(Key::ForwardStrandReadDepths.to_string(), "ADF");
        assert_eq!(Key::ReverseStrandReadDepths.to_string(), "ADR");
        assert_eq!(Key::AlleleFrequencies.to_string(), "AF");
        assert_eq!(Key::TotalAlleleCount.to_string(), "AN");
        assert_eq!(Key::BaseQuality.to_string(), "BQ");
        assert_eq!(Key::Cigar.to_string(), "CIGAR");
        assert_eq!(Key::IsInDbSnp.to_string(), "DB");
        assert_eq!(Key::TotalDepth.to_string(), "DP");
        // assert_eq!(Key::EndPosition.to_string(), "END");
        assert_eq!(Key::IsInHapMap2.to_string(), "H2");
        assert_eq!(Key::IsInHapMap3.to_string(), "H3");
        assert_eq!(Key::MappingQuality.to_string(), "MQ");
        assert_eq!(Key::ZeroMappingQualityCount.to_string(), "MQ0");
        assert_eq!(Key::SamplesWithDataCount.to_string(), "NS");
        assert_eq!(Key::StrandBias.to_string(), "SB");
        assert_eq!(Key::IsSomaticMutation.to_string(), "SOMATIC");
        assert_eq!(Key::IsValidated.to_string(), "VALIDATED");
        assert_eq!(Key::IsIn1000Genomes.to_string(), "1000G");

        assert_eq!(Key::IsImprecise.to_string(), "IMPRECISE");
        assert_eq!(Key::IsNovel.to_string(), "NOVEL");
        assert_eq!(Key::EndPosition.to_string(), "END");
        assert_eq!(Key::SvType.to_string(), "SVTYPE");
        assert_eq!(Key::SvLengths.to_string(), "SVLEN");
        assert_eq!(Key::PositionConfidenceIntervals.to_string(), "CIPOS");
        assert_eq!(Key::EndConfidenceIntervals.to_string(), "CIEND");
        assert_eq!(Key::MicrohomologyLengths.to_string(), "HOMLEN");
        assert_eq!(Key::MicrohomologySequences.to_string(), "HOMSEQ");
        assert_eq!(Key::BreakpointIds.to_string(), "BKPTID");
        assert_eq!(Key::MobileElementInfo.to_string(), "MEINFO");
        assert_eq!(Key::MobileElementTransductionInfo.to_string(), "METRANS");
        assert_eq!(Key::DbvId.to_string(), "DGVID");
        assert_eq!(Key::DbVarId.to_string(), "DBVARID");
        assert_eq!(Key::DbRipId.to_string(), "DBRIPID");
        assert_eq!(Key::MateBreakendIds.to_string(), "MATEID");
        assert_eq!(Key::PartnerBreakendId.to_string(), "PARID");
        assert_eq!(Key::BreakendEventId.to_string(), "EVENT");
        assert_eq!(Key::BreakendConfidenceIntervals.to_string(), "CILEN");
        // assert_eq!(Key::BreakendReadDepth.to_string(), "DP");
        assert_eq!(Key::AdjacentReadDepths.to_string(), "DPADJ");
        assert_eq!(Key::BreakendCopyNumber.to_string(), "CN");
        assert_eq!(Key::AdjacentCopyNumber.to_string(), "CNADJ");
        assert_eq!(Key::CopyNumberConfidenceIntervals.to_string(), "CICN");
        assert_eq!(
            Key::AdjacentCopyNumberConfidenceIntervals.to_string(),
            "CICNADJ"
        );

        assert_eq!(Key::Other(String::from("NDLS")).to_string(), "NDLS");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("AA".parse(), Ok(Key::AncestralAllele));
        assert_eq!("AC".parse(), Ok(Key::AlleleCount));
        assert_eq!("AD".parse(), Ok(Key::TotalReadDepths));
        assert_eq!("ADF".parse(), Ok(Key::ForwardStrandReadDepths));
        assert_eq!("ADR".parse(), Ok(Key::ReverseStrandReadDepths));
        assert_eq!("AF".parse(), Ok(Key::AlleleFrequencies));
        assert_eq!("AN".parse(), Ok(Key::TotalAlleleCount));
        assert_eq!("BQ".parse(), Ok(Key::BaseQuality));
        assert_eq!("CIGAR".parse(), Ok(Key::Cigar));
        assert_eq!("DB".parse(), Ok(Key::IsInDbSnp));
        assert_eq!("DP".parse(), Ok(Key::TotalDepth));
        // assert_eq!("END".parse(), Ok(Key::EndPosition));
        assert_eq!("H2".parse(), Ok(Key::IsInHapMap2));
        assert_eq!("H3".parse(), Ok(Key::IsInHapMap3));
        assert_eq!("MQ".parse(), Ok(Key::MappingQuality));
        assert_eq!("MQ0".parse(), Ok(Key::ZeroMappingQualityCount));
        assert_eq!("NS".parse(), Ok(Key::SamplesWithDataCount));
        assert_eq!("SB".parse(), Ok(Key::StrandBias));
        assert_eq!("SOMATIC".parse(), Ok(Key::IsSomaticMutation));
        assert_eq!("VALIDATED".parse(), Ok(Key::IsValidated));
        assert_eq!("1000G".parse(), Ok(Key::IsIn1000Genomes));

        assert_eq!("IMPRECISE".parse(), Ok(Key::IsImprecise));
        assert_eq!("NOVEL".parse(), Ok(Key::IsNovel));
        assert_eq!("END".parse(), Ok(Key::EndPosition));
        assert_eq!("SVTYPE".parse(), Ok(Key::SvType));
        assert_eq!("SVLEN".parse(), Ok(Key::SvLengths));
        assert_eq!("CIPOS".parse(), Ok(Key::PositionConfidenceIntervals));
        assert_eq!("CIEND".parse(), Ok(Key::EndConfidenceIntervals));
        assert_eq!("HOMLEN".parse(), Ok(Key::MicrohomologyLengths));
        assert_eq!("HOMSEQ".parse(), Ok(Key::MicrohomologySequences));
        assert_eq!("BKPTID".parse(), Ok(Key::BreakpointIds));
        assert_eq!("MEINFO".parse(), Ok(Key::MobileElementInfo));
        assert_eq!("METRANS".parse(), Ok(Key::MobileElementTransductionInfo));
        assert_eq!("DGVID".parse(), Ok(Key::DbvId));
        assert_eq!("DBVARID".parse(), Ok(Key::DbVarId));
        assert_eq!("DBRIPID".parse(), Ok(Key::DbRipId));
        assert_eq!("MATEID".parse(), Ok(Key::MateBreakendIds));
        assert_eq!("PARID".parse(), Ok(Key::PartnerBreakendId));
        assert_eq!("EVENT".parse(), Ok(Key::BreakendEventId));
        assert_eq!("CILEN".parse(), Ok(Key::BreakendConfidenceIntervals));
        // assert_eq!("DP".parse(), Ok(Key::BreakendReadDepth));
        assert_eq!("DPADJ".parse(), Ok(Key::AdjacentReadDepths));
        assert_eq!("CN".parse(), Ok(Key::BreakendCopyNumber));
        assert_eq!("CNADJ".parse(), Ok(Key::AdjacentCopyNumber));
        assert_eq!("CICN".parse(), Ok(Key::CopyNumberConfidenceIntervals));
        assert_eq!(
            "CICNADJ".parse(),
            Ok(Key::AdjacentCopyNumberConfidenceIntervals)
        );

        assert_eq!("NDLS".parse(), Ok(Key::Other(String::from("NDLS"))));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
        assert_eq!("8D".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!(".N".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!("A!".parse::<Key>(), Err(ParseError::Invalid));
    }

    #[test]
    fn test_number() {
        assert_eq!(number(&Key::AncestralAllele), Some(Number::Count(1)));
        assert_eq!(number(&Key::AlleleCount), Some(Number::A));
        assert_eq!(number(&Key::TotalReadDepths), Some(Number::R));
        assert_eq!(number(&Key::ForwardStrandReadDepths), Some(Number::R));
        assert_eq!(number(&Key::ReverseStrandReadDepths), Some(Number::R));
        assert_eq!(number(&Key::AlleleFrequencies), Some(Number::A));
        assert_eq!(number(&Key::TotalAlleleCount), Some(Number::Count(1)));
        assert_eq!(number(&Key::BaseQuality), Some(Number::Count(1)));
        assert_eq!(number(&Key::Cigar), Some(Number::A));
        assert_eq!(number(&Key::IsInDbSnp), Some(Number::Count(0)));
        assert_eq!(number(&Key::TotalDepth), Some(Number::Count(1)));
        // assert_eq!(number(&Key::EndPosition), Some(Number::Count(1)));
        assert_eq!(number(&Key::IsInHapMap2), Some(Number::Count(0)));
        assert_eq!(number(&Key::IsInHapMap3), Some(Number::Count(0)));
        assert_eq!(number(&Key::MappingQuality), Some(Number::Count(1)));
        assert_eq!(
            number(&Key::ZeroMappingQualityCount),
            Some(Number::Count(1))
        );
        assert_eq!(number(&Key::SamplesWithDataCount), Some(Number::Count(1)));
        assert_eq!(number(&Key::StrandBias), Some(Number::Count(4)));
        assert_eq!(number(&Key::IsSomaticMutation), Some(Number::Count(0)));
        assert_eq!(number(&Key::IsValidated), Some(Number::Count(0)));
        assert_eq!(number(&Key::IsIn1000Genomes), Some(Number::Count(0)));

        assert_eq!(number(&Key::IsImprecise), Some(Number::Count(0)));
        assert_eq!(number(&Key::IsNovel), Some(Number::Count(0)));
        assert_eq!(number(&Key::EndPosition), Some(Number::Count(1)));
        assert_eq!(number(&Key::SvType), Some(Number::Count(1)));
        assert_eq!(number(&Key::SvLengths), Some(Number::Unknown));
        assert_eq!(
            number(&Key::PositionConfidenceIntervals),
            Some(Number::Count(2))
        );
        assert_eq!(number(&Key::EndConfidenceIntervals), Some(Number::Count(2)));
        assert_eq!(number(&Key::MicrohomologyLengths), Some(Number::Unknown));
        assert_eq!(number(&Key::MicrohomologySequences), Some(Number::Unknown));
        assert_eq!(number(&Key::BreakpointIds), Some(Number::Unknown));
        assert_eq!(number(&Key::MobileElementInfo), Some(Number::Count(4)));
        assert_eq!(
            number(&Key::MobileElementTransductionInfo),
            Some(Number::Count(4))
        );
        assert_eq!(number(&Key::DbvId), Some(Number::Count(1)));
        assert_eq!(number(&Key::DbVarId), Some(Number::Count(1)));
        assert_eq!(number(&Key::DbRipId), Some(Number::Count(1)));
        assert_eq!(number(&Key::MateBreakendIds), Some(Number::Unknown));
        assert_eq!(number(&Key::PartnerBreakendId), Some(Number::Count(1)));
        assert_eq!(number(&Key::BreakendEventId), Some(Number::Count(1)));
        assert_eq!(
            number(&Key::BreakendConfidenceIntervals),
            Some(Number::Count(2))
        );
        // assert_eq!(number(&Key::BreakendReadDepth), Some(Number::Count(1)));
        assert_eq!(number(&Key::AdjacentReadDepths), Some(Number::Unknown));
        assert_eq!(number(&Key::BreakendCopyNumber), Some(Number::Count(1)));
        assert_eq!(number(&Key::AdjacentCopyNumber), Some(Number::Unknown));
        assert_eq!(
            number(&Key::CopyNumberConfidenceIntervals),
            Some(Number::Count(2))
        );
        assert_eq!(
            number(&Key::AdjacentCopyNumberConfidenceIntervals),
            Some(Number::Unknown)
        );

        assert!(number(&Key::Other(String::from("NDLS"))).is_none());
    }

    #[test]
    fn test_ty() {
        assert_eq!(ty(&Key::AncestralAllele), Some(Type::String));
        assert_eq!(ty(&Key::AlleleCount), Some(Type::Integer));
        assert_eq!(ty(&Key::TotalReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::ForwardStrandReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::ReverseStrandReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::AlleleFrequencies), Some(Type::Float));
        assert_eq!(ty(&Key::TotalAlleleCount), Some(Type::Integer));
        assert_eq!(ty(&Key::BaseQuality), Some(Type::Float));
        assert_eq!(ty(&Key::Cigar), Some(Type::String));
        assert_eq!(ty(&Key::IsInDbSnp), Some(Type::Flag));
        assert_eq!(ty(&Key::TotalDepth), Some(Type::Integer));
        // assert_eq!(ty(&Key::EndPosition), Some(Type::Integer));
        assert_eq!(ty(&Key::IsInHapMap2), Some(Type::Flag));
        assert_eq!(ty(&Key::IsInHapMap3), Some(Type::Flag));
        assert_eq!(ty(&Key::MappingQuality), Some(Type::Float));
        assert_eq!(ty(&Key::ZeroMappingQualityCount), Some(Type::Integer));
        assert_eq!(ty(&Key::SamplesWithDataCount), Some(Type::Integer));
        assert_eq!(ty(&Key::StrandBias), Some(Type::Integer));
        assert_eq!(ty(&Key::IsSomaticMutation), Some(Type::Flag));
        assert_eq!(ty(&Key::IsValidated), Some(Type::Flag));
        assert_eq!(ty(&Key::IsIn1000Genomes), Some(Type::Flag));

        assert_eq!(ty(&Key::IsImprecise), Some(Type::Flag));
        assert_eq!(ty(&Key::IsNovel), Some(Type::Flag));
        assert_eq!(ty(&Key::EndPosition), Some(Type::Integer));
        assert_eq!(ty(&Key::SvType), Some(Type::String));
        assert_eq!(ty(&Key::SvLengths), Some(Type::Integer));
        assert_eq!(ty(&Key::PositionConfidenceIntervals), Some(Type::Integer));
        assert_eq!(ty(&Key::EndConfidenceIntervals), Some(Type::Integer));
        assert_eq!(ty(&Key::MicrohomologyLengths), Some(Type::Integer));
        assert_eq!(ty(&Key::MicrohomologySequences), Some(Type::String));
        assert_eq!(ty(&Key::BreakpointIds), Some(Type::String));
        assert_eq!(ty(&Key::MobileElementInfo), Some(Type::String));
        assert_eq!(ty(&Key::MobileElementTransductionInfo), Some(Type::String));
        assert_eq!(ty(&Key::DbvId), Some(Type::String));
        assert_eq!(ty(&Key::DbVarId), Some(Type::String));
        assert_eq!(ty(&Key::DbRipId), Some(Type::String));
        assert_eq!(ty(&Key::MateBreakendIds), Some(Type::String));
        assert_eq!(ty(&Key::PartnerBreakendId), Some(Type::String));
        assert_eq!(ty(&Key::BreakendEventId), Some(Type::String));
        assert_eq!(ty(&Key::BreakendConfidenceIntervals), Some(Type::Integer));
        // assert_eq!(ty(&Key::BreakendReadDepth), Some(Type::Integer));
        assert_eq!(ty(&Key::AdjacentReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::BreakendCopyNumber), Some(Type::Integer));
        assert_eq!(ty(&Key::AdjacentCopyNumber), Some(Type::Integer));
        assert_eq!(ty(&Key::CopyNumberConfidenceIntervals), Some(Type::Integer));
        assert_eq!(
            ty(&Key::AdjacentCopyNumberConfidenceIntervals),
            Some(Type::Integer)
        );

        assert!(ty(&Key::Other(String::from("NDLS"))).is_none());
    }

    #[test]
    fn test_description() {
        assert_eq!(description(&Key::AncestralAllele), Some("Ancestral allele"));
        assert_eq!(
            description(&Key::AlleleCount),
            Some("Allele count in genotypes, for each ALT allele, in the same order as listed")
        );
        assert_eq!(
            description(&Key::TotalReadDepths),
            Some("Total read depth for each allele")
        );
        assert_eq!(
            description(&Key::ForwardStrandReadDepths),
            Some("Read depth for each allele on the forward strand")
        );
        assert_eq!(
            description(&Key::ReverseStrandReadDepths),
            Some("Read depth for each allele on the reverse strand")
        );
        assert_eq!(
            description(&Key::AlleleFrequencies),
            Some("Allele frequency for each ALT allele in the same order as listed")
        );
        assert_eq!(
            description(&Key::TotalAlleleCount),
            Some("Total number of alleles in called genotypes")
        );
        assert_eq!(description(&Key::BaseQuality), Some("RMS base quality"));
        assert_eq!(
            description(&Key::Cigar),
            Some(
                "Cigar string describing how to align an alternate allele to the reference allele"
            )
        );
        assert_eq!(description(&Key::IsInDbSnp), Some("dbSNP membership"));
        assert_eq!(
            description(&Key::TotalDepth),
            Some("Combined depth across samples")
        );
        // Self::EndPosition.description(), Some("End position on CHROM"));
        assert_eq!(description(&Key::IsInHapMap2), Some("HapMap2 membership"));
        assert_eq!(description(&Key::IsInHapMap3), Some("HapMap3 membership"));
        assert_eq!(
            description(&Key::MappingQuality),
            Some("RMS mapping quality")
        );
        assert_eq!(
            description(&Key::ZeroMappingQualityCount),
            Some("Number of MAPQ == 0 reads")
        );
        assert_eq!(
            description(&Key::SamplesWithDataCount),
            Some("Number of samples with data")
        );
        assert_eq!(description(&Key::StrandBias), Some("Strand bias"));
        assert_eq!(
            description(&Key::IsSomaticMutation),
            Some("Somatic mutation")
        );
        assert_eq!(
            description(&Key::IsValidated),
            Some("Validated by follow-up experiment")
        );
        assert_eq!(
            description(&Key::IsIn1000Genomes),
            Some("1000 Genomes membership")
        );

        assert_eq!(
            description(&Key::IsImprecise),
            Some("Imprecise structural variation")
        );
        assert_eq!(
            description(&Key::IsNovel),
            Some("Indicates a novel structural variation")
        );
        assert_eq!(
            description(&Key::EndPosition),
            Some("End position of the variant described in this record")
        );
        assert_eq!(
            description(&Key::SvType),
            Some("Type of structural variant")
        );
        assert_eq!(
            description(&Key::SvLengths),
            Some("Difference in length between REF and ALT alleles")
        );
        assert_eq!(
            description(&Key::PositionConfidenceIntervals),
            Some("Confidence interval around POS for imprecise variants")
        );
        assert_eq!(
            description(&Key::EndConfidenceIntervals),
            Some("Confidence interval around END for imprecise variants")
        );
        assert_eq!(
            description(&Key::MicrohomologyLengths),
            Some("Length of base pair identical micro-homology at event breakpoints")
        );
        assert_eq!(
            description(&Key::MicrohomologySequences),
            Some("Sequence of base pair identical micro-homology at event breakpoints")
        );
        assert_eq!(
            description(&Key::BreakpointIds),
            Some("ID of the assembled alternate allele in the assembly file")
        );
        assert_eq!(
            description(&Key::MobileElementInfo),
            Some("Mobile element info of the form NAME,START,END,POLARITY")
        );
        assert_eq!(
            description(&Key::MobileElementTransductionInfo),
            Some("Mobile element transduction info of the form CHR,START,END,POLARITY")
        );
        assert_eq!(
            description(&Key::DbvId),
            Some("ID of this element in Database of Genomic Variation")
        );
        assert_eq!(
            description(&Key::DbVarId),
            Some("ID of this element in DBVAR")
        );
        assert_eq!(
            description(&Key::DbRipId),
            Some("ID of this element in DBRIP")
        );
        assert_eq!(
            description(&Key::MateBreakendIds),
            Some("ID of mate breakends")
        );
        assert_eq!(
            description(&Key::PartnerBreakendId),
            Some("ID of partner breakend")
        );
        assert_eq!(
            description(&Key::BreakendEventId),
            Some("ID of event associated to breakend")
        );
        assert_eq!(
            description(&Key::BreakendConfidenceIntervals),
            Some("Confidence interval around the inserted material between breakends")
        );
        // Self::BreakendReadDepth.description(), Some("Read Depth of segment containing breakend"));
        assert_eq!(
            description(&Key::AdjacentReadDepths),
            Some("Read Depth of adjacency")
        );
        assert_eq!(
            description(&Key::BreakendCopyNumber),
            Some("Copy number of segment containing breakend")
        );
        assert_eq!(
            description(&Key::AdjacentCopyNumber),
            Some("Copy number of adjacency")
        );
        assert_eq!(
            description(&Key::CopyNumberConfidenceIntervals),
            Some("Confidence interval around copy number for the segment")
        );
        assert_eq!(
            description(&Key::AdjacentCopyNumberConfidenceIntervals),
            Some("Confidence interval around copy number for the adjacency")
        );

        assert!(description(&Key::Other(String::from("NDLS"))).is_none());
    }
}
