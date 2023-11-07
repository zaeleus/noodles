//! VCF record info key.

mod other;
mod standard;

use std::{
    borrow::Borrow,
    error, fmt,
    hash::{Hash, Hasher},
    str::FromStr,
};

pub use self::{other::Other, standard::Standard};
use crate::header::FileFormat;

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
///
/// Deprecated in VCF 4.4.
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

/// Type of associated event (`EVENTTYPE`).
///
/// Added in VCF 4.4.
pub const EVENT_TYPE: Key = Key::Standard(Standard::EventType);

/// Confidence interval around the inserted material between breakends (`CILEN`).
pub const BREAKEND_CONFIDENCE_INTERVALS: Key = Key::Standard(Standard::BreakendConfidenceIntervals);

/// Read Depth of adjacency (`DPADJ`).
///
/// Removed in VCF 4.4.
pub const ADJACENT_READ_DEPTHS: Key = Key::Standard(Standard::AdjacentReadDepths);

/// Copy number of segment containing breakend (`CN`).
pub const BREAKEND_COPY_NUMBER: Key = Key::Standard(Standard::BreakendCopyNumber);

/// Copy number of adjacency (`CNADJ`).
///
/// Removed in VCF 4.4.
pub const ADJACENT_COPY_NUMBER: Key = Key::Standard(Standard::AdjacentCopyNumber);

/// Confidence interval around copy number for the segment (`CICN`).
pub const COPY_NUMBER_CONFIDENCE_INTERVALS: Key =
    Key::Standard(Standard::CopyNumberConfidenceIntervals);

/// Confidence interval around copy number for the adjacency (`CICNADJ`).
///
/// Removed in VCF 4.4.
pub const ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS: Key =
    Key::Standard(Standard::AdjacentCopyNumberConfidenceIntervals);

/// Copy number of adjacency (`CNADJ`).
///
/// Added in VCF 4.4.
pub const SV_CLAIM: Key = Key::Standard(Standard::SvClaim);

/// Total number of repeat sequences in this allele (`RN`).
///
/// Added in VCF 4.4.
pub const TOTAL_REPEAT_SEQUENCE_COUNTS: Key = Key::Standard(Standard::TotalRepeatSequenceCounts);

/// Repeat unit sequence of the corresponding repeat sequence (`RUS`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_SEQUENCES: Key = Key::Standard(Standard::RepeatUnitSequences);

/// Repeat unit length of the corresponding repeat sequence (`RUL`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_LENGTHS: Key = Key::Standard(Standard::RepeatUnitLengths);

/// Repeat unit count of corresponding repeat sequence (`RUC`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_COUNTS: Key = Key::Standard(Standard::RepeatUnitCounts);

/// Total number of bases in the corresponding repeat sequence (`RB`).
///
/// Added in VCF 4.4.
pub const TOTAL_REPEAT_SEQUENCE_BASE_COUNTS: Key =
    Key::Standard(Standard::TotalRepeatSequenceBaseCounts);

/// Confidence interval around RUC (`CIRUC`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_COUNT_CONFIDENCE_INTERVALS: Key =
    Key::Standard(Standard::RepeatUnitCountConfidenceIntervals);

/// Confidence interval around RB (`CIRB`).
///
/// Added in VCF 4.4.
pub const TOTAL_REPEAT_SEQUENCE_BASE_COUNT_CONFIDENCE_INTERVALS: Key =
    Key::Standard(Standard::TotalRepeatSequenceBaseCountConfidenceIntervals);

/// Number of bases in each individual repeat unit (`RUB`).
///
/// Added in VCF 4.4.
pub const REPEAT_UNIT_BASE_COUNTS: Key = Key::Standard(Standard::RepeatUnitBaseCounts);

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

/// A VCF header info key.
#[derive(Clone, Debug, Eq, PartialEq)]
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

impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_ref().hash(state);
    }
}

impl TryFrom<(FileFormat, &str)> for Key {
    type Error = ParseError;

    fn try_from((file_format, s): (FileFormat, &str)) -> Result<Self, Self::Error> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        s.parse()
            .map(Self::Standard)
            .or_else(|_| Other::try_from((file_format, s)).map(Self::Other))
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
        assert_eq!(SV_CLAIM.to_string(), "SVCLAIM");
        assert_eq!(TOTAL_REPEAT_SEQUENCE_COUNTS.to_string(), "RN");
        assert_eq!(REPEAT_UNIT_SEQUENCES.to_string(), "RUS");
        assert_eq!(REPEAT_UNIT_LENGTHS.to_string(), "RUL");
        assert_eq!(REPEAT_UNIT_COUNTS.to_string(), "RUC");
        assert_eq!(TOTAL_REPEAT_SEQUENCE_BASE_COUNTS.to_string(), "RB");
        assert_eq!(REPEAT_UNIT_COUNT_CONFIDENCE_INTERVALS.to_string(), "CIRUC");
        assert_eq!(
            TOTAL_REPEAT_SEQUENCE_BASE_COUNT_CONFIDENCE_INTERVALS.to_string(),
            "CIRB"
        );
        assert_eq!(REPEAT_UNIT_BASE_COUNTS.to_string(), "RUB");

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
        assert_eq!("EVENTTYPE".parse(), Ok(EVENT_TYPE));
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
        assert_eq!("SVCLAIM".parse(), Ok(SV_CLAIM));
        assert_eq!("RN".parse(), Ok(TOTAL_REPEAT_SEQUENCE_COUNTS));
        assert_eq!("RUS".parse(), Ok(REPEAT_UNIT_SEQUENCES));
        assert_eq!("RUL".parse(), Ok(REPEAT_UNIT_LENGTHS));
        assert_eq!("RUC".parse(), Ok(REPEAT_UNIT_COUNTS));
        assert_eq!("RB".parse(), Ok(TOTAL_REPEAT_SEQUENCE_BASE_COUNTS));
        assert_eq!("CIRUC".parse(), Ok(REPEAT_UNIT_COUNT_CONFIDENCE_INTERVALS));
        assert_eq!(
            "CIRB".parse(),
            Ok(TOTAL_REPEAT_SEQUENCE_BASE_COUNT_CONFIDENCE_INTERVALS)
        );
        assert_eq!("RUB".parse(), Ok(REPEAT_UNIT_BASE_COUNTS));

        assert_eq!("NDLS".parse(), Ok(Key::Other(Other(String::from("NDLS")))));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
        assert_eq!("8D".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!(".N".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!("A!".parse::<Key>(), Err(ParseError::Invalid));
    }
}
