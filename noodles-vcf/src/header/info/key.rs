//! VCF header info key.

mod v4_3;
mod v4_4;

use crate::header::{record::value::map::info::Type, Number};

use std::{borrow::Borrow, error, fmt, str::FromStr};

/// Ancestral allele (`AA`).
pub const ANCESTRAL_ALLELE: Key = Key::Standard(Standard::AncestralAllele);

/// Allele count in genotypes, for each ALT allele, in the same order as listed (`AC`).
pub const ALLELE_COUNT: Key = Key::Standard(Standard::AlleleCount);

/// Total read depth for each allele (`AD`).
pub const TOTAL_READ_DEPTHS: Key = Key::Standard(Standard::TotalReadDepths);

/// Read depth for each allele on the forward strand (`ADF`).
pub const FORWARD_STRAND_READ_DEPTHS: Key = Key::Standard(Standard::ForwardStrandReadDepths);

/// Read depth for each allele on the reverse strand (`ADR`).
pub const REVERSE_STRAND_READ_DEPTHS: Key = Key::Standard(Standard::ReverseStrandReadDepths);

/// Allele frequency for each ALT allele in the same order as listed (`AF`).
pub const ALLELE_FREQUENCIES: Key = Key::Standard(Standard::AlleleFrequencies);

/// Total number of alleles in called genotypes (`AN`).
pub const TOTAL_ALLELE_COUNT: Key = Key::Standard(Standard::TotalAlleleCount);

/// RMS base quality (`BQ`).
pub const BASE_QUALITY: Key = Key::Standard(Standard::BaseQuality);

/// Cigar string describing how to align an alternate allele to the reference allele (`CIGAR`).
pub const CIGAR: Key = Key::Standard(Standard::Cigar);

/// dbSNP membership (`DB`).
pub const IS_IN_DB_SNP: Key = Key::Standard(Standard::IsInDbSnp);

/// Combined depth across samples (`DP`).
pub const TOTAL_DEPTH: Key = Key::Standard(Standard::TotalDepth);

/// HapMap2 membership (`H2`).
pub const IS_IN_HAP_MAP_2: Key = Key::Standard(Standard::IsInHapMap2);

/// HapMap3 membership (`H3`).
pub const IS_IN_HAP_MAP_3: Key = Key::Standard(Standard::IsInHapMap3);

/// RMS mapping quality (`MQ`).
pub const MAPPING_QUALITY: Key = Key::Standard(Standard::MappingQuality);

/// Number of MAPQ == 0 reads (`MQ0`).
pub const ZERO_MAPPING_QUALITY_COUNT: Key = Key::Standard(Standard::ZeroMappingQualityCount);

/// Number of samples with data (`NS`).
pub const SAMPLES_WITH_DATA_COUNT: Key = Key::Standard(Standard::SamplesWithDataCount);

/// Strand bias (`SB`).
pub const STRAND_BIAS: Key = Key::Standard(Standard::StrandBias);

/// Somatic mutation (`SOMATIC`).
pub const IS_SOMATIC_MUTATION: Key = Key::Standard(Standard::IsSomaticMutation);

/// Validated by follow-up experiment (`VALIDATED`).
pub const IS_VALIDATED: Key = Key::Standard(Standard::IsValidated);

/// 1000 Genomes membership (`1000G`).
pub const IS_IN_1000_GENOMES: Key = Key::Standard(Standard::IsIn1000Genomes);

/// Imprecise structural variation (`IMPRECISE`).
pub const IS_IMPRECISE: Key = Key::Standard(Standard::IsImprecise);

/// Indicates a novel structural variation (`NOVEL`).
pub const IS_NOVEL: Key = Key::Standard(Standard::IsNovel);

/// End position of the variant described in this record (`END`).
pub const END_POSITION: Key = Key::Standard(Standard::EndPosition);

/// Type of structural variant (`SVTYPE`).
pub const SV_TYPE: Key = Key::Standard(Standard::SvType);

/// Difference in length between REF and ALT alleles (`SVLEN`).
pub const SV_LENGTHS: Key = Key::Standard(Standard::SvLengths);

/// Confidence interval around POS for imprecise variants (`CIPOS`).
pub const POSITION_CONFIDENCE_INTERVALS: Key = Key::Standard(Standard::PositionConfidenceIntervals);

/// Confidence interval around END for imprecise variants (`CIEND`).
pub const END_CONFIDENCE_INTERVALS: Key = Key::Standard(Standard::EndConfidenceIntervals);

/// Length of base pair identical micro-homology at event breakpoints (`HOMLEN`).
pub const MICROHOMOLOGY_LENGTHS: Key = Key::Standard(Standard::MicrohomologyLengths);

/// Sequence of base pair identical micro-homology at event breakpoints (`HOMSEQ`).
pub const MICROHOMOLOGY_SEQUENCES: Key = Key::Standard(Standard::MicrohomologySequences);

/// ID of the assembled alternate allele in the assembly file (`BKPTID`).
pub const BREAKPOINT_IDS: Key = Key::Standard(Standard::BreakpointIds);

/// Mobile element info of the form NAME,START,END,POLARITY (`MEINFO`).
pub const MOBILE_ELEMENT_INFO: Key = Key::Standard(Standard::MobileElementInfo);

/// Mobile element transduction info of the form CHR,START,END,POLARITY (`METRANS`).
pub const MOBILE_ELEMENT_TRANSDUCTION_INFO: Key =
    Key::Standard(Standard::MobileElementTransductionInfo);

/// ID of this element in Database of Genomic Variation (`DBVID`).
pub const DBV_ID: Key = Key::Standard(Standard::DbvId);

/// ID of this element in DBVAR (`DBVARID`).
pub const DB_VAR_ID: Key = Key::Standard(Standard::DbVarId);

/// ID of this element in DBRIP (`DBRIPID`).
pub const DB_RIP_ID: Key = Key::Standard(Standard::DbRipId);

/// ID of mate breakends (`MATEID`).
pub const MATE_BREAKEND_IDS: Key = Key::Standard(Standard::MateBreakendIds);

/// ID of partner breakend (`PARID`).
pub const PARTNER_BREAKEND_ID: Key = Key::Standard(Standard::PartnerBreakendId);

/// ID of event associated to breakend (`EVENT`).
pub const BREAKEND_EVENT_ID: Key = Key::Standard(Standard::BreakendEventId);

/// Confidence interval around the inserted material between breakends (`CILEN`).
pub const BREAKEND_CONFIDENCE_INTERVALS: Key = Key::Standard(Standard::BreakendConfidenceIntervals);

/// Read Depth of adjacency (`DPADJ`).
pub const ADJACENT_READ_DEPTHS: Key = Key::Standard(Standard::AdjacentReadDepths);

/// Copy number of segment containing breakend (`CN`).
pub const BREAKEND_COPY_NUMBER: Key = Key::Standard(Standard::BreakendCopyNumber);

/// Copy number of adjacency (`CNADJ`).
pub const ADJACENT_COPY_NUMBER: Key = Key::Standard(Standard::AdjacentCopyNumber);

/// Confidence interval around copy number for the segment (`CICN`).
pub const COPY_NUMBER_CONFIDENCE_INTERVALS: Key =
    Key::Standard(Standard::CopyNumberConfidenceIntervals);

/// Confidence interval around copy number for the adjacency (`CICNADJ`).
pub const ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS: Key =
    Key::Standard(Standard::AdjacentCopyNumberConfidenceIntervals);

/// A VCF header info key.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Standard {
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
            Self::BreakendConfidenceIntervals => "CILEN",
            // Self::BreakendReadDepth => "DP",
            Self::AdjacentReadDepths => "DPADJ",
            Self::BreakendCopyNumber => "CN",
            Self::AdjacentCopyNumber => "CNADJ",
            Self::CopyNumberConfidenceIntervals => "CICN",
            Self::AdjacentCopyNumberConfidenceIntervals => "CICNADJ",
        }
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
            "CILEN" => Ok(Self::BreakendConfidenceIntervals),
            // "DP" => Ok(Self::BreakendReadDepth),
            "DPADJ" => Ok(Self::AdjacentReadDepths),
            "CN" => Ok(Self::BreakendCopyNumber),
            "CNADJ" => Ok(Self::AdjacentCopyNumber),
            "CICN" => Ok(Self::CopyNumberConfidenceIntervals),
            "CICNADJ" => Ok(Self::AdjacentCopyNumberConfidenceIntervals),

            _ => Err(ParseError::Invalid),
        }
    }
}

/// A non-reserved VCF header info key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Other(String);

impl AsRef<str> for Other {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl FromStr for Other {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if is_valid_name(s) {
            Ok(Self(s.into()))
        } else {
            Err(ParseError::Invalid)
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

/// A VCF header info key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Key {
    /// A reserved key.
    Standard(Standard),
    /// Any other non-reserved key.
    Other(Other),
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::Standard(k) => k.as_ref(),
            Self::Other(k) => k.as_ref(),
        }
    }
}

impl Borrow<str> for Key {
    fn borrow(&self) -> &str {
        self.as_ref()
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        s.parse()
            .map(Self::Standard)
            .or_else(|_| s.parse().map(Self::Other))
    }
}

pub(crate) fn number(key: &Key) -> Option<Number> {
    match key {
        Key::Standard(k) => v4_3::definition(*k).map(|(number, _, _)| number),
        Key::Other(_) => None,
    }
}

pub(crate) fn ty(key: &Key) -> Option<Type> {
    match key {
        Key::Standard(k) => v4_3::definition(*k).map(|(_, ty, _)| ty),
        Key::Other(_) => None,
    }
}

pub(crate) fn description(key: &Key) -> Option<&str> {
    match key {
        Key::Standard(k) => v4_3::definition(*k).map(|(_, _, description)| description),
        Key::Other(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(ANCESTRAL_ALLELE.to_string(), "AA");
        assert_eq!(ALLELE_COUNT.to_string(), "AC");
        assert_eq!(TOTAL_READ_DEPTHS.to_string(), "AD");
        assert_eq!(FORWARD_STRAND_READ_DEPTHS.to_string(), "ADF");
        assert_eq!(REVERSE_STRAND_READ_DEPTHS.to_string(), "ADR");
        assert_eq!(ALLELE_FREQUENCIES.to_string(), "AF");
        assert_eq!(TOTAL_ALLELE_COUNT.to_string(), "AN");
        assert_eq!(BASE_QUALITY.to_string(), "BQ");
        assert_eq!(CIGAR.to_string(), "CIGAR");
        assert_eq!(IS_IN_DB_SNP.to_string(), "DB");
        assert_eq!(TOTAL_DEPTH.to_string(), "DP");
        // assert_eq!(END_POSITION.to_string(), "END");
        assert_eq!(IS_IN_HAP_MAP_2.to_string(), "H2");
        assert_eq!(IS_IN_HAP_MAP_3.to_string(), "H3");
        assert_eq!(MAPPING_QUALITY.to_string(), "MQ");
        assert_eq!(ZERO_MAPPING_QUALITY_COUNT.to_string(), "MQ0");
        assert_eq!(SAMPLES_WITH_DATA_COUNT.to_string(), "NS");
        assert_eq!(STRAND_BIAS.to_string(), "SB");
        assert_eq!(IS_SOMATIC_MUTATION.to_string(), "SOMATIC");
        assert_eq!(IS_VALIDATED.to_string(), "VALIDATED");
        assert_eq!(IS_IN_1000_GENOMES.to_string(), "1000G");

        assert_eq!(IS_IMPRECISE.to_string(), "IMPRECISE");
        assert_eq!(IS_NOVEL.to_string(), "NOVEL");
        assert_eq!(END_POSITION.to_string(), "END");
        assert_eq!(SV_TYPE.to_string(), "SVTYPE");
        assert_eq!(SV_LENGTHS.to_string(), "SVLEN");
        assert_eq!(POSITION_CONFIDENCE_INTERVALS.to_string(), "CIPOS");
        assert_eq!(END_CONFIDENCE_INTERVALS.to_string(), "CIEND");
        assert_eq!(MICROHOMOLOGY_LENGTHS.to_string(), "HOMLEN");
        assert_eq!(MICROHOMOLOGY_SEQUENCES.to_string(), "HOMSEQ");
        assert_eq!(BREAKPOINT_IDS.to_string(), "BKPTID");
        assert_eq!(MOBILE_ELEMENT_INFO.to_string(), "MEINFO");
        assert_eq!(MOBILE_ELEMENT_TRANSDUCTION_INFO.to_string(), "METRANS");
        assert_eq!(DBV_ID.to_string(), "DGVID");
        assert_eq!(DB_VAR_ID.to_string(), "DBVARID");
        assert_eq!(DB_RIP_ID.to_string(), "DBRIPID");
        assert_eq!(MATE_BREAKEND_IDS.to_string(), "MATEID");
        assert_eq!(PARTNER_BREAKEND_ID.to_string(), "PARID");
        assert_eq!(BREAKEND_EVENT_ID.to_string(), "EVENT");
        assert_eq!(BREAKEND_CONFIDENCE_INTERVALS.to_string(), "CILEN");
        // assert_eq!(Key::BreakendReadDepth.to_string(), "DP");
        assert_eq!(ADJACENT_READ_DEPTHS.to_string(), "DPADJ");
        assert_eq!(BREAKEND_COPY_NUMBER.to_string(), "CN");
        assert_eq!(ADJACENT_COPY_NUMBER.to_string(), "CNADJ");
        assert_eq!(COPY_NUMBER_CONFIDENCE_INTERVALS.to_string(), "CICN");
        assert_eq!(
            ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS.to_string(),
            "CICNADJ"
        );

        assert_eq!(Key::Other(Other(String::from("NDLS"))).to_string(), "NDLS");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("AA".parse(), Ok(ANCESTRAL_ALLELE));
        assert_eq!("AC".parse(), Ok(ALLELE_COUNT));
        assert_eq!("AD".parse(), Ok(TOTAL_READ_DEPTHS));
        assert_eq!("ADF".parse(), Ok(FORWARD_STRAND_READ_DEPTHS));
        assert_eq!("ADR".parse(), Ok(REVERSE_STRAND_READ_DEPTHS));
        assert_eq!("AF".parse(), Ok(ALLELE_FREQUENCIES));
        assert_eq!("AN".parse(), Ok(TOTAL_ALLELE_COUNT));
        assert_eq!("BQ".parse(), Ok(BASE_QUALITY));
        assert_eq!("CIGAR".parse(), Ok(CIGAR));
        assert_eq!("DB".parse(), Ok(IS_IN_DB_SNP));
        assert_eq!("DP".parse(), Ok(TOTAL_DEPTH));
        // assert_eq!("END".parse(), Ok(END_POSITION));
        assert_eq!("H2".parse(), Ok(IS_IN_HAP_MAP_2));
        assert_eq!("H3".parse(), Ok(IS_IN_HAP_MAP_3));
        assert_eq!("MQ".parse(), Ok(MAPPING_QUALITY));
        assert_eq!("MQ0".parse(), Ok(ZERO_MAPPING_QUALITY_COUNT));
        assert_eq!("NS".parse(), Ok(SAMPLES_WITH_DATA_COUNT));
        assert_eq!("SB".parse(), Ok(STRAND_BIAS));
        assert_eq!("SOMATIC".parse(), Ok(IS_SOMATIC_MUTATION));
        assert_eq!("VALIDATED".parse(), Ok(IS_VALIDATED));
        assert_eq!("1000G".parse(), Ok(IS_IN_1000_GENOMES));

        assert_eq!("IMPRECISE".parse(), Ok(IS_IMPRECISE));
        assert_eq!("NOVEL".parse(), Ok(IS_NOVEL));
        assert_eq!("END".parse(), Ok(END_POSITION));
        assert_eq!("SVTYPE".parse(), Ok(SV_TYPE));
        assert_eq!("SVLEN".parse(), Ok(SV_LENGTHS));
        assert_eq!("CIPOS".parse(), Ok(POSITION_CONFIDENCE_INTERVALS));
        assert_eq!("CIEND".parse(), Ok(END_CONFIDENCE_INTERVALS));
        assert_eq!("HOMLEN".parse(), Ok(MICROHOMOLOGY_LENGTHS));
        assert_eq!("HOMSEQ".parse(), Ok(MICROHOMOLOGY_SEQUENCES));
        assert_eq!("BKPTID".parse(), Ok(BREAKPOINT_IDS));
        assert_eq!("MEINFO".parse(), Ok(MOBILE_ELEMENT_INFO));
        assert_eq!("METRANS".parse(), Ok(MOBILE_ELEMENT_TRANSDUCTION_INFO));
        assert_eq!("DGVID".parse(), Ok(DBV_ID));
        assert_eq!("DBVARID".parse(), Ok(DB_VAR_ID));
        assert_eq!("DBRIPID".parse(), Ok(DB_RIP_ID));
        assert_eq!("MATEID".parse(), Ok(MATE_BREAKEND_IDS));
        assert_eq!("PARID".parse(), Ok(PARTNER_BREAKEND_ID));
        assert_eq!("EVENT".parse(), Ok(BREAKEND_EVENT_ID));
        assert_eq!("CILEN".parse(), Ok(BREAKEND_CONFIDENCE_INTERVALS));
        // assert_eq!("DP".parse(), Ok(Key::BreakendReadDepth));
        assert_eq!("DPADJ".parse(), Ok(ADJACENT_READ_DEPTHS));
        assert_eq!("CN".parse(), Ok(BREAKEND_COPY_NUMBER));
        assert_eq!("CNADJ".parse(), Ok(ADJACENT_COPY_NUMBER));
        assert_eq!("CICN".parse(), Ok(COPY_NUMBER_CONFIDENCE_INTERVALS));
        assert_eq!(
            "CICNADJ".parse(),
            Ok(ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS)
        );

        assert_eq!("NDLS".parse(), Ok(Key::Other(Other(String::from("NDLS")))));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
        assert_eq!("8D".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!(".N".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!("A!".parse::<Key>(), Err(ParseError::Invalid));
    }

    #[test]
    fn test_number() {
        assert_eq!(number(&ANCESTRAL_ALLELE), Some(Number::Count(1)));
        assert_eq!(number(&ALLELE_COUNT), Some(Number::A));
        assert_eq!(number(&TOTAL_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&FORWARD_STRAND_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&REVERSE_STRAND_READ_DEPTHS), Some(Number::R));
        assert_eq!(number(&ALLELE_FREQUENCIES), Some(Number::A));
        assert_eq!(number(&TOTAL_ALLELE_COUNT), Some(Number::Count(1)));
        assert_eq!(number(&BASE_QUALITY), Some(Number::Count(1)));
        assert_eq!(number(&CIGAR), Some(Number::A));
        assert_eq!(number(&IS_IN_DB_SNP), Some(Number::Count(0)));
        assert_eq!(number(&TOTAL_DEPTH), Some(Number::Count(1)));
        // assert_eq!(number(&END_POSITION), Some(Number::Count(1)));
        assert_eq!(number(&IS_IN_HAP_MAP_2), Some(Number::Count(0)));
        assert_eq!(number(&IS_IN_HAP_MAP_3), Some(Number::Count(0)));
        assert_eq!(number(&MAPPING_QUALITY), Some(Number::Count(1)));
        assert_eq!(number(&ZERO_MAPPING_QUALITY_COUNT), Some(Number::Count(1)));
        assert_eq!(number(&SAMPLES_WITH_DATA_COUNT), Some(Number::Count(1)));
        assert_eq!(number(&STRAND_BIAS), Some(Number::Count(4)));
        assert_eq!(number(&IS_SOMATIC_MUTATION), Some(Number::Count(0)));
        assert_eq!(number(&IS_VALIDATED), Some(Number::Count(0)));
        assert_eq!(number(&IS_IN_1000_GENOMES), Some(Number::Count(0)));

        assert_eq!(number(&IS_IMPRECISE), Some(Number::Count(0)));
        assert_eq!(number(&IS_NOVEL), Some(Number::Count(0)));
        assert_eq!(number(&END_POSITION), Some(Number::Count(1)));
        assert_eq!(number(&SV_TYPE), Some(Number::Count(1)));
        assert_eq!(number(&SV_LENGTHS), Some(Number::Unknown));
        assert_eq!(
            number(&POSITION_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        assert_eq!(number(&END_CONFIDENCE_INTERVALS), Some(Number::Count(2)));
        assert_eq!(number(&MICROHOMOLOGY_LENGTHS), Some(Number::Unknown));
        assert_eq!(number(&MICROHOMOLOGY_SEQUENCES), Some(Number::Unknown));
        assert_eq!(number(&BREAKPOINT_IDS), Some(Number::Unknown));
        assert_eq!(number(&MOBILE_ELEMENT_INFO), Some(Number::Count(4)));
        assert_eq!(
            number(&MOBILE_ELEMENT_TRANSDUCTION_INFO),
            Some(Number::Count(4))
        );
        assert_eq!(number(&DBV_ID), Some(Number::Count(1)));
        assert_eq!(number(&DB_VAR_ID), Some(Number::Count(1)));
        assert_eq!(number(&DB_RIP_ID), Some(Number::Count(1)));
        assert_eq!(number(&MATE_BREAKEND_IDS), Some(Number::Unknown));
        assert_eq!(number(&PARTNER_BREAKEND_ID), Some(Number::Count(1)));
        assert_eq!(number(&BREAKEND_EVENT_ID), Some(Number::Count(1)));
        assert_eq!(
            number(&BREAKEND_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        // assert_eq!(number(&Key::BreakendReadDepth), Some(Number::Count(1)));
        assert_eq!(number(&ADJACENT_READ_DEPTHS), Some(Number::Unknown));
        assert_eq!(number(&BREAKEND_COPY_NUMBER), Some(Number::Count(1)));
        assert_eq!(number(&ADJACENT_COPY_NUMBER), Some(Number::Unknown));
        assert_eq!(
            number(&COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Number::Count(2))
        );
        assert_eq!(
            number(&ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Number::Unknown)
        );

        assert!(number(&Key::Other(Other(String::from("NDLS")))).is_none());
    }

    #[test]
    fn test_ty() {
        assert_eq!(ty(&ANCESTRAL_ALLELE), Some(Type::String));
        assert_eq!(ty(&ALLELE_COUNT), Some(Type::Integer));
        assert_eq!(ty(&TOTAL_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&FORWARD_STRAND_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&REVERSE_STRAND_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&ALLELE_FREQUENCIES), Some(Type::Float));
        assert_eq!(ty(&TOTAL_ALLELE_COUNT), Some(Type::Integer));
        assert_eq!(ty(&BASE_QUALITY), Some(Type::Float));
        assert_eq!(ty(&CIGAR), Some(Type::String));
        assert_eq!(ty(&IS_IN_DB_SNP), Some(Type::Flag));
        assert_eq!(ty(&TOTAL_DEPTH), Some(Type::Integer));
        // assert_eq!(ty(&END_POSITION), Some(Type::Integer));
        assert_eq!(ty(&IS_IN_HAP_MAP_2), Some(Type::Flag));
        assert_eq!(ty(&IS_IN_HAP_MAP_3), Some(Type::Flag));
        assert_eq!(ty(&MAPPING_QUALITY), Some(Type::Float));
        assert_eq!(ty(&ZERO_MAPPING_QUALITY_COUNT), Some(Type::Integer));
        assert_eq!(ty(&SAMPLES_WITH_DATA_COUNT), Some(Type::Integer));
        assert_eq!(ty(&STRAND_BIAS), Some(Type::Integer));
        assert_eq!(ty(&IS_SOMATIC_MUTATION), Some(Type::Flag));
        assert_eq!(ty(&IS_VALIDATED), Some(Type::Flag));
        assert_eq!(ty(&IS_IN_1000_GENOMES), Some(Type::Flag));

        assert_eq!(ty(&IS_IMPRECISE), Some(Type::Flag));
        assert_eq!(ty(&IS_NOVEL), Some(Type::Flag));
        assert_eq!(ty(&END_POSITION), Some(Type::Integer));
        assert_eq!(ty(&SV_TYPE), Some(Type::String));
        assert_eq!(ty(&SV_LENGTHS), Some(Type::Integer));
        assert_eq!(ty(&POSITION_CONFIDENCE_INTERVALS), Some(Type::Integer));
        assert_eq!(ty(&END_CONFIDENCE_INTERVALS), Some(Type::Integer));
        assert_eq!(ty(&MICROHOMOLOGY_LENGTHS), Some(Type::Integer));
        assert_eq!(ty(&MICROHOMOLOGY_SEQUENCES), Some(Type::String));
        assert_eq!(ty(&BREAKPOINT_IDS), Some(Type::String));
        assert_eq!(ty(&MOBILE_ELEMENT_INFO), Some(Type::String));
        assert_eq!(ty(&MOBILE_ELEMENT_TRANSDUCTION_INFO), Some(Type::String));
        assert_eq!(ty(&DBV_ID), Some(Type::String));
        assert_eq!(ty(&DB_VAR_ID), Some(Type::String));
        assert_eq!(ty(&DB_RIP_ID), Some(Type::String));
        assert_eq!(ty(&MATE_BREAKEND_IDS), Some(Type::String));
        assert_eq!(ty(&PARTNER_BREAKEND_ID), Some(Type::String));
        assert_eq!(ty(&BREAKEND_EVENT_ID), Some(Type::String));
        assert_eq!(ty(&BREAKEND_CONFIDENCE_INTERVALS), Some(Type::Integer));
        // assert_eq!(ty(&Key::BreakendReadDepth), Some(Type::Integer));
        assert_eq!(ty(&ADJACENT_READ_DEPTHS), Some(Type::Integer));
        assert_eq!(ty(&BREAKEND_COPY_NUMBER), Some(Type::Integer));
        assert_eq!(ty(&ADJACENT_COPY_NUMBER), Some(Type::Integer));
        assert_eq!(ty(&COPY_NUMBER_CONFIDENCE_INTERVALS), Some(Type::Integer));
        assert_eq!(
            ty(&ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some(Type::Integer)
        );

        assert!(ty(&Key::Other(Other(String::from("NDLS")))).is_none());
    }

    #[test]
    fn test_description() {
        assert_eq!(description(&ANCESTRAL_ALLELE), Some("Ancestral allele"));
        assert_eq!(
            description(&ALLELE_COUNT),
            Some("Allele count in genotypes, for each ALT allele, in the same order as listed")
        );
        assert_eq!(
            description(&TOTAL_READ_DEPTHS),
            Some("Total read depth for each allele")
        );
        assert_eq!(
            description(&FORWARD_STRAND_READ_DEPTHS),
            Some("Read depth for each allele on the forward strand")
        );
        assert_eq!(
            description(&REVERSE_STRAND_READ_DEPTHS),
            Some("Read depth for each allele on the reverse strand")
        );
        assert_eq!(
            description(&ALLELE_FREQUENCIES),
            Some("Allele frequency for each ALT allele in the same order as listed")
        );
        assert_eq!(
            description(&TOTAL_ALLELE_COUNT),
            Some("Total number of alleles in called genotypes")
        );
        assert_eq!(description(&BASE_QUALITY), Some("RMS base quality"));
        assert_eq!(
            description(&CIGAR),
            Some(
                "Cigar string describing how to align an alternate allele to the reference allele"
            )
        );
        assert_eq!(description(&IS_IN_DB_SNP), Some("dbSNP membership"));
        assert_eq!(
            description(&TOTAL_DEPTH),
            Some("Combined depth across samples")
        );
        // Self::EndPosition.description(), Some("End position on CHROM"));
        assert_eq!(description(&IS_IN_HAP_MAP_2), Some("HapMap2 membership"));
        assert_eq!(description(&IS_IN_HAP_MAP_3), Some("HapMap3 membership"));
        assert_eq!(description(&MAPPING_QUALITY), Some("RMS mapping quality"));
        assert_eq!(
            description(&ZERO_MAPPING_QUALITY_COUNT),
            Some("Number of MAPQ == 0 reads")
        );
        assert_eq!(
            description(&SAMPLES_WITH_DATA_COUNT),
            Some("Number of samples with data")
        );
        assert_eq!(description(&STRAND_BIAS), Some("Strand bias"));
        assert_eq!(description(&IS_SOMATIC_MUTATION), Some("Somatic mutation"));
        assert_eq!(
            description(&IS_VALIDATED),
            Some("Validated by follow-up experiment")
        );
        assert_eq!(
            description(&IS_IN_1000_GENOMES),
            Some("1000 Genomes membership")
        );

        assert_eq!(
            description(&IS_IMPRECISE),
            Some("Imprecise structural variation")
        );
        assert_eq!(
            description(&IS_NOVEL),
            Some("Indicates a novel structural variation")
        );
        assert_eq!(
            description(&END_POSITION),
            Some("End position of the variant described in this record")
        );
        assert_eq!(description(&SV_TYPE), Some("Type of structural variant"));
        assert_eq!(
            description(&SV_LENGTHS),
            Some("Difference in length between REF and ALT alleles")
        );
        assert_eq!(
            description(&POSITION_CONFIDENCE_INTERVALS),
            Some("Confidence interval around POS for imprecise variants")
        );
        assert_eq!(
            description(&END_CONFIDENCE_INTERVALS),
            Some("Confidence interval around END for imprecise variants")
        );
        assert_eq!(
            description(&MICROHOMOLOGY_LENGTHS),
            Some("Length of base pair identical micro-homology at event breakpoints")
        );
        assert_eq!(
            description(&MICROHOMOLOGY_SEQUENCES),
            Some("Sequence of base pair identical micro-homology at event breakpoints")
        );
        assert_eq!(
            description(&BREAKPOINT_IDS),
            Some("ID of the assembled alternate allele in the assembly file")
        );
        assert_eq!(
            description(&MOBILE_ELEMENT_INFO),
            Some("Mobile element info of the form NAME,START,END,POLARITY")
        );
        assert_eq!(
            description(&MOBILE_ELEMENT_TRANSDUCTION_INFO),
            Some("Mobile element transduction info of the form CHR,START,END,POLARITY")
        );
        assert_eq!(
            description(&DBV_ID),
            Some("ID of this element in Database of Genomic Variation")
        );
        assert_eq!(description(&DB_VAR_ID), Some("ID of this element in DBVAR"));
        assert_eq!(description(&DB_RIP_ID), Some("ID of this element in DBRIP"));
        assert_eq!(
            description(&MATE_BREAKEND_IDS),
            Some("ID of mate breakends")
        );
        assert_eq!(
            description(&PARTNER_BREAKEND_ID),
            Some("ID of partner breakend")
        );
        assert_eq!(
            description(&BREAKEND_EVENT_ID),
            Some("ID of event associated to breakend")
        );
        assert_eq!(
            description(&BREAKEND_CONFIDENCE_INTERVALS),
            Some("Confidence interval around the inserted material between breakends")
        );
        // Self::BreakendReadDepth.description(), Some("Read Depth of segment containing breakend"));
        assert_eq!(
            description(&ADJACENT_READ_DEPTHS),
            Some("Read Depth of adjacency")
        );
        assert_eq!(
            description(&BREAKEND_COPY_NUMBER),
            Some("Copy number of segment containing breakend")
        );
        assert_eq!(
            description(&ADJACENT_COPY_NUMBER),
            Some("Copy number of adjacency")
        );
        assert_eq!(
            description(&COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some("Confidence interval around copy number for the segment")
        );
        assert_eq!(
            description(&ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS),
            Some("Confidence interval around copy number for the adjacency")
        );

        assert!(description(&Key::Other(Other(String::from("NDLS")))).is_none());
    }
}
