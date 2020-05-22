use crate::header::{info::Type, Number};

use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Key {
    // § 1.6.1 Fixed Fields (2020-04-02)
    AncestralAllele,
    AlleleCount,
    TotalReadDepths,
    ForwardStrandReadDepths,
    ReverseStrandReadDepths,
    AlleleFrequencies,
    TotalAlleleCount,
    BaseQuality,
    Cigar,
    IsInDbSnp,
    TotalDepth,
    // EndPosition,
    IsInHapMap2,
    IsInHapMap3,
    MappingQuality,
    ZeroMappingQualityCount,
    SamplesWithDataCount,
    StrandBias,
    IsSomaticMutation,
    IsValidated,
    IsIn1000Genomes,

    // § 3 INFO keys used for structural variants (2020-04-02)
    IsImprecise,
    IsNovel,
    EndPosition,
    SvType,
    SvLengths,
    PositionConfidenceIntervals,
    EndConfidenceIntervals,
    MicrohomologyLengths,
    MicrohomologySequences,
    BreakpointIds,
    MobileElementInfo,
    MobileElementTransductionInfo,
    DbvId,
    DbVarId,
    DbRipId,
    MateBreakendIds,
    PartnerBreakendId,
    BreakendEventId,
    BreakendConfidenceIntervals,
    // BreakendReadDepth,
    AdjacentReadDepths,
    BreakendCopyNumber,
    AdjacentCopyNumber,
    CopyNumberConfidenceIntervals,
    AdjacentCopyNumberConfidenceIntervals,

    Other(String, Number, Type),
}

impl Key {
    pub fn number(&self) -> Number {
        match self {
            Self::AncestralAllele => Number::Count(1),
            Self::AlleleCount => Number::A,
            Self::TotalReadDepths => Number::R,
            Self::ForwardStrandReadDepths => Number::R,
            Self::ReverseStrandReadDepths => Number::R,
            Self::AlleleFrequencies => Number::A,
            Self::TotalAlleleCount => Number::Count(1),
            Self::BaseQuality => Number::Count(1),
            Self::Cigar => Number::A,
            Self::IsInDbSnp => Number::Count(0),
            Self::TotalDepth => Number::Count(1),
            // Self::EndPosition => Number::Count(1),
            Self::IsInHapMap2 => Number::Count(0),
            Self::IsInHapMap3 => Number::Count(0),
            Self::MappingQuality => Number::Count(1),
            Self::ZeroMappingQualityCount => Number::Count(1),
            Self::SamplesWithDataCount => Number::Count(1),
            Self::StrandBias => Number::Count(4),
            Self::IsSomaticMutation => Number::Count(0),
            Self::IsValidated => Number::Count(0),
            Self::IsIn1000Genomes => Number::Count(0),

            Self::IsImprecise => Number::Count(0),
            Self::IsNovel => Number::Count(0),
            Self::EndPosition => Number::Count(1),
            Self::SvType => Number::Count(1),
            Self::SvLengths => Number::Unknown,
            Self::PositionConfidenceIntervals => Number::Count(2),
            Self::EndConfidenceIntervals => Number::Count(2),
            Self::MicrohomologyLengths => Number::Unknown,
            Self::MicrohomologySequences => Number::Unknown,
            Self::BreakpointIds => Number::Unknown,
            Self::MobileElementInfo => Number::Count(4),
            Self::MobileElementTransductionInfo => Number::Count(4),
            Self::DbvId => Number::Count(1),
            Self::DbVarId => Number::Count(1),
            Self::DbRipId => Number::Count(1),
            Self::MateBreakendIds => Number::Unknown,
            Self::PartnerBreakendId => Number::Count(1),
            Self::BreakendEventId => Number::Count(1),
            Self::BreakendConfidenceIntervals => Number::Count(2),
            // Self::BreakendReadDepth => Number::Count(1),
            Self::AdjacentReadDepths => Number::Unknown,
            Self::BreakendCopyNumber => Number::Count(1),
            Self::AdjacentCopyNumber => Number::Unknown,
            Self::CopyNumberConfidenceIntervals => Number::Count(2),
            Self::AdjacentCopyNumberConfidenceIntervals => Number::Unknown,

            Self::Other(_, number, _) => *number,
        }
    }

    pub fn ty(&self) -> Type {
        match self {
            Self::AncestralAllele => Type::String,
            Self::AlleleCount => Type::Integer,
            Self::TotalReadDepths => Type::Integer,
            Self::ForwardStrandReadDepths => Type::Integer,
            Self::ReverseStrandReadDepths => Type::Integer,
            Self::AlleleFrequencies => Type::Float,
            Self::TotalAlleleCount => Type::Integer,
            Self::BaseQuality => Type::Float,
            Self::Cigar => Type::String,
            Self::IsInDbSnp => Type::Flag,
            Self::TotalDepth => Type::Integer,
            // Self::EndPosition => Type::Integer,
            Self::IsInHapMap2 => Type::Flag,
            Self::IsInHapMap3 => Type::Flag,
            Self::MappingQuality => Type::Float,
            Self::ZeroMappingQualityCount => Type::Integer,
            Self::SamplesWithDataCount => Type::Integer,
            Self::StrandBias => Type::Integer,
            Self::IsSomaticMutation => Type::Flag,
            Self::IsValidated => Type::Flag,
            Self::IsIn1000Genomes => Type::Flag,

            Self::IsImprecise => Type::Flag,
            Self::IsNovel => Type::Flag,
            Self::EndPosition => Type::Integer,
            Self::SvType => Type::String,
            Self::SvLengths => Type::Integer,
            Self::PositionConfidenceIntervals => Type::Integer,
            Self::EndConfidenceIntervals => Type::Integer,
            Self::MicrohomologyLengths => Type::Integer,
            Self::MicrohomologySequences => Type::String,
            Self::BreakpointIds => Type::String,
            Self::MobileElementInfo => Type::String,
            Self::MobileElementTransductionInfo => Type::String,
            Self::DbvId => Type::String,
            Self::DbVarId => Type::String,
            Self::DbRipId => Type::String,
            Self::MateBreakendIds => Type::String,
            Self::PartnerBreakendId => Type::String,
            Self::BreakendEventId => Type::String,
            Self::BreakendConfidenceIntervals => Type::Integer,
            // Self::BreakendReadDepth => Type::Integer,
            Self::AdjacentReadDepths => Type::Integer,
            Self::BreakendCopyNumber => Type::Integer,
            Self::AdjacentCopyNumber => Type::Integer,
            Self::CopyNumberConfidenceIntervals => Type::Integer,
            Self::AdjacentCopyNumberConfidenceIntervals => Type::Integer,

            Self::Other(_, _, ty) => *ty,
        }
    }
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

            Self::Other(key, ..) => key,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid info key: {}", self.0)
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError(s.into()));
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

            _ => Ok(Self::Other(s.into(), Number::Count(1), Type::String)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_number() {
        assert_eq!(Key::AncestralAllele.number(), Number::Count(1));
        assert_eq!(Key::AlleleCount.number(), Number::A);
        assert_eq!(Key::TotalReadDepths.number(), Number::R);
        assert_eq!(Key::ForwardStrandReadDepths.number(), Number::R);
        assert_eq!(Key::ReverseStrandReadDepths.number(), Number::R);
        assert_eq!(Key::AlleleFrequencies.number(), Number::A);
        assert_eq!(Key::TotalAlleleCount.number(), Number::Count(1));
        assert_eq!(Key::BaseQuality.number(), Number::Count(1));
        assert_eq!(Key::Cigar.number(), Number::A);
        assert_eq!(Key::IsInDbSnp.number(), Number::Count(0));
        assert_eq!(Key::TotalDepth.number(), Number::Count(1));
        // assert_eq!(Key::EndPosition.number(), Number::Count(1));
        assert_eq!(Key::IsInHapMap2.number(), Number::Count(0));
        assert_eq!(Key::IsInHapMap3.number(), Number::Count(0));
        assert_eq!(Key::MappingQuality.number(), Number::Count(1));
        assert_eq!(Key::ZeroMappingQualityCount.number(), Number::Count(1));
        assert_eq!(Key::SamplesWithDataCount.number(), Number::Count(1));
        assert_eq!(Key::StrandBias.number(), Number::Count(4));
        assert_eq!(Key::IsSomaticMutation.number(), Number::Count(0));
        assert_eq!(Key::IsValidated.number(), Number::Count(0));
        assert_eq!(Key::IsIn1000Genomes.number(), Number::Count(0));

        assert_eq!(Key::IsImprecise.number(), Number::Count(0));
        assert_eq!(Key::IsNovel.number(), Number::Count(0));
        assert_eq!(Key::EndPosition.number(), Number::Count(1));
        assert_eq!(Key::SvType.number(), Number::Count(1));
        assert_eq!(Key::SvLengths.number(), Number::Unknown);
        assert_eq!(Key::PositionConfidenceIntervals.number(), Number::Count(2));
        assert_eq!(Key::EndConfidenceIntervals.number(), Number::Count(2));
        assert_eq!(Key::MicrohomologyLengths.number(), Number::Unknown);
        assert_eq!(Key::MicrohomologySequences.number(), Number::Unknown);
        assert_eq!(Key::BreakpointIds.number(), Number::Unknown);
        assert_eq!(Key::MobileElementInfo.number(), Number::Count(4));
        assert_eq!(
            Key::MobileElementTransductionInfo.number(),
            Number::Count(4)
        );
        assert_eq!(Key::DbvId.number(), Number::Count(1));
        assert_eq!(Key::DbVarId.number(), Number::Count(1));
        assert_eq!(Key::DbRipId.number(), Number::Count(1));
        assert_eq!(Key::MateBreakendIds.number(), Number::Unknown);
        assert_eq!(Key::PartnerBreakendId.number(), Number::Count(1));
        assert_eq!(Key::BreakendEventId.number(), Number::Count(1));
        assert_eq!(Key::BreakendConfidenceIntervals.number(), Number::Count(2));
        // assert_eq!(Key::BreakendReadDepth.number(), Number::Count(1));
        assert_eq!(Key::AdjacentReadDepths.number(), Number::Unknown);
        assert_eq!(Key::BreakendCopyNumber.number(), Number::Count(1));
        assert_eq!(Key::AdjacentCopyNumber.number(), Number::Unknown);
        assert_eq!(
            Key::CopyNumberConfidenceIntervals.number(),
            Number::Count(2)
        );
        assert_eq!(
            Key::AdjacentCopyNumberConfidenceIntervals.number(),
            Number::Unknown
        );

        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).number(),
            Number::Count(1)
        );
    }

    #[test]
    fn test_ty() {
        assert_eq!(Key::AncestralAllele.ty(), Type::String);
        assert_eq!(Key::AlleleCount.ty(), Type::Integer);
        assert_eq!(Key::TotalReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ForwardStrandReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ReverseStrandReadDepths.ty(), Type::Integer);
        assert_eq!(Key::AlleleFrequencies.ty(), Type::Float);
        assert_eq!(Key::TotalAlleleCount.ty(), Type::Integer);
        assert_eq!(Key::BaseQuality.ty(), Type::Float);
        assert_eq!(Key::Cigar.ty(), Type::String);
        assert_eq!(Key::IsInDbSnp.ty(), Type::Flag);
        assert_eq!(Key::TotalDepth.ty(), Type::Integer);
        // assert_eq!(Key::EndPosition.ty(), Type::Integer);
        assert_eq!(Key::IsInHapMap2.ty(), Type::Flag);
        assert_eq!(Key::IsInHapMap3.ty(), Type::Flag);
        assert_eq!(Key::MappingQuality.ty(), Type::Float);
        assert_eq!(Key::ZeroMappingQualityCount.ty(), Type::Integer);
        assert_eq!(Key::SamplesWithDataCount.ty(), Type::Integer);
        assert_eq!(Key::StrandBias.ty(), Type::Integer);
        assert_eq!(Key::IsSomaticMutation.ty(), Type::Flag);
        assert_eq!(Key::IsValidated.ty(), Type::Flag);
        assert_eq!(Key::IsIn1000Genomes.ty(), Type::Flag);

        assert_eq!(Key::IsImprecise.ty(), Type::Flag);
        assert_eq!(Key::IsNovel.ty(), Type::Flag);
        assert_eq!(Key::EndPosition.ty(), Type::Integer);
        assert_eq!(Key::SvType.ty(), Type::String);
        assert_eq!(Key::SvLengths.ty(), Type::Integer);
        assert_eq!(Key::PositionConfidenceIntervals.ty(), Type::Integer);
        assert_eq!(Key::EndConfidenceIntervals.ty(), Type::Integer);
        assert_eq!(Key::MicrohomologyLengths.ty(), Type::Integer);
        assert_eq!(Key::MicrohomologySequences.ty(), Type::String);
        assert_eq!(Key::BreakpointIds.ty(), Type::String);
        assert_eq!(Key::MobileElementInfo.ty(), Type::String);
        assert_eq!(Key::MobileElementTransductionInfo.ty(), Type::String);
        assert_eq!(Key::DbvId.ty(), Type::String);
        assert_eq!(Key::DbVarId.ty(), Type::String);
        assert_eq!(Key::DbRipId.ty(), Type::String);
        assert_eq!(Key::MateBreakendIds.ty(), Type::String);
        assert_eq!(Key::PartnerBreakendId.ty(), Type::String);
        assert_eq!(Key::BreakendEventId.ty(), Type::String);
        assert_eq!(Key::BreakendConfidenceIntervals.ty(), Type::Integer);
        // assert_eq!(Key::BreakendReadDepth.ty(), Type::Integer);
        assert_eq!(Key::AdjacentReadDepths.ty(), Type::Integer);
        assert_eq!(Key::BreakendCopyNumber.ty(), Type::Integer);
        assert_eq!(Key::AdjacentCopyNumber.ty(), Type::Integer);
        assert_eq!(Key::CopyNumberConfidenceIntervals.ty(), Type::Integer);
        assert_eq!(
            Key::AdjacentCopyNumberConfidenceIntervals.ty(),
            Type::Integer
        );

        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).ty(),
            Type::String
        );
    }

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

        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).to_string(),
            "NDLS"
        );
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("AA".parse::<Key>()?, Key::AncestralAllele);
        assert_eq!("AC".parse::<Key>()?, Key::AlleleCount);
        assert_eq!("AD".parse::<Key>()?, Key::TotalReadDepths);
        assert_eq!("ADF".parse::<Key>()?, Key::ForwardStrandReadDepths);
        assert_eq!("ADR".parse::<Key>()?, Key::ReverseStrandReadDepths);
        assert_eq!("AF".parse::<Key>()?, Key::AlleleFrequencies);
        assert_eq!("AN".parse::<Key>()?, Key::TotalAlleleCount);
        assert_eq!("BQ".parse::<Key>()?, Key::BaseQuality);
        assert_eq!("CIGAR".parse::<Key>()?, Key::Cigar);
        assert_eq!("DB".parse::<Key>()?, Key::IsInDbSnp);
        assert_eq!("DP".parse::<Key>()?, Key::TotalDepth);
        // assert_eq!("END".parse::<Key>()?, Key::EndPosition);
        assert_eq!("H2".parse::<Key>()?, Key::IsInHapMap2);
        assert_eq!("H3".parse::<Key>()?, Key::IsInHapMap3);
        assert_eq!("MQ".parse::<Key>()?, Key::MappingQuality);
        assert_eq!("MQ0".parse::<Key>()?, Key::ZeroMappingQualityCount);
        assert_eq!("NS".parse::<Key>()?, Key::SamplesWithDataCount);
        assert_eq!("SB".parse::<Key>()?, Key::StrandBias);
        assert_eq!("SOMATIC".parse::<Key>()?, Key::IsSomaticMutation);
        assert_eq!("VALIDATED".parse::<Key>()?, Key::IsValidated);
        assert_eq!("1000G".parse::<Key>()?, Key::IsIn1000Genomes);

        assert_eq!("IMPRECISE".parse::<Key>()?, Key::IsImprecise);
        assert_eq!("NOVEL".parse::<Key>()?, Key::IsNovel);
        assert_eq!("END".parse::<Key>()?, Key::EndPosition);
        assert_eq!("SVTYPE".parse::<Key>()?, Key::SvType);
        assert_eq!("SVLEN".parse::<Key>()?, Key::SvLengths);
        assert_eq!("CIPOS".parse::<Key>()?, Key::PositionConfidenceIntervals);
        assert_eq!("CIEND".parse::<Key>()?, Key::EndConfidenceIntervals);
        assert_eq!("HOMLEN".parse::<Key>()?, Key::MicrohomologyLengths);
        assert_eq!("HOMSEQ".parse::<Key>()?, Key::MicrohomologySequences);
        assert_eq!("BKPTID".parse::<Key>()?, Key::BreakpointIds);
        assert_eq!("MEINFO".parse::<Key>()?, Key::MobileElementInfo);
        assert_eq!(
            "METRANS".parse::<Key>()?,
            Key::MobileElementTransductionInfo
        );
        assert_eq!("DGVID".parse::<Key>()?, Key::DbvId);
        assert_eq!("DBVARID".parse::<Key>()?, Key::DbVarId);
        assert_eq!("DBRIPID".parse::<Key>()?, Key::DbRipId);
        assert_eq!("MATEID".parse::<Key>()?, Key::MateBreakendIds);
        assert_eq!("PARID".parse::<Key>()?, Key::PartnerBreakendId);
        assert_eq!("EVENT".parse::<Key>()?, Key::BreakendEventId);
        assert_eq!("CILEN".parse::<Key>()?, Key::BreakendConfidenceIntervals);
        // assert_eq!("DP".parse::<Key>()?, Key::BreakendReadDepth);
        assert_eq!("DPADJ".parse::<Key>()?, Key::AdjacentReadDepths);
        assert_eq!("CN".parse::<Key>()?, Key::BreakendCopyNumber);
        assert_eq!("CNADJ".parse::<Key>()?, Key::AdjacentCopyNumber);
        assert_eq!("CICN".parse::<Key>()?, Key::CopyNumberConfidenceIntervals);
        assert_eq!(
            "CICNADJ".parse::<Key>()?,
            Key::AdjacentCopyNumberConfidenceIntervals
        );

        assert_eq!(
            "NDLS".parse::<Key>()?,
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String)
        );

        assert!("".parse::<Key>().is_err());

        Ok(())
    }
}
