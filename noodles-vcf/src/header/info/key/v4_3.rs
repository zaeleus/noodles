use crate::{
    header::{record::value::map::info::Type, Number},
    record::info::field::key::Standard,
};

pub(super) fn definition(key: Standard) -> Option<(Number, Type, &'static str)> {
    match key {
        Standard::AncestralAllele => Some((Number::Count(1), Type::String, "Ancestral allele")),
        Standard::AlleleCount => Some((
            Number::A,
            Type::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        )),
        Standard::TotalReadDepths => {
            Some((Number::R, Type::Integer, "Total read depth for each allele"))
        }
        Standard::ForwardStrandReadDepths => Some((
            Number::R,
            Type::Integer,
            "Read depth for each allele on the forward strand",
        )),
        Standard::ReverseStrandReadDepths => Some((
            Number::R,
            Type::Integer,
            "Read depth for each allele on the reverse strand",
        )),
        Standard::AlleleFrequencies => Some((
            Number::A,
            Type::Float,
            "Allele frequency for each ALT allele in the same order as listed",
        )),
        Standard::TotalAlleleCount => Some((
            Number::Count(1),
            Type::Integer,
            "Total number of alleles in called genotypes",
        )),
        Standard::BaseQuality => Some((Number::Count(1), Type::Float, "RMS base quality")),
        Standard::Cigar => Some((
            Number::A,
            Type::String,
            "Cigar string describing how to align an alternate allele to the reference allele",
        )),
        Standard::IsInDbSnp => Some((Number::Count(0), Type::Flag, "dbSNP membership")),
        Standard::TotalDepth => Some((
            Number::Count(1),
            Type::Integer,
            "Combined depth across samples",
        )),
        Standard::IsInHapMap2 => Some((Number::Count(0), Type::Flag, "HapMap2 membership")),
        Standard::IsInHapMap3 => Some((Number::Count(0), Type::Flag, "HapMap3 membership")),
        Standard::MappingQuality => Some((Number::Count(1), Type::Float, "RMS mapping quality")),
        Standard::ZeroMappingQualityCount => {
            Some((Number::Count(1), Type::Integer, "Number of MAPQ == 0 reads"))
        }
        Standard::SamplesWithDataCount => Some((
            Number::Count(1),
            Type::Integer,
            "Number of samples with data",
        )),
        Standard::StrandBias => Some((Number::Count(4), Type::Integer, "Strand bias")),
        Standard::IsSomaticMutation => Some((Number::Count(0), Type::Flag, "Somatic mutation")),
        Standard::IsValidated => Some((
            Number::Count(0),
            Type::Flag,
            "Validated by follow-up experiment",
        )),
        Standard::IsIn1000Genomes => {
            Some((Number::Count(0), Type::Flag, "1000 Genomes membership"))
        }

        Standard::IsImprecise => Some((
            Number::Count(0),
            Type::Flag,
            "Imprecise structural variation",
        )),
        Standard::IsNovel => Some((
            Number::Count(0),
            Type::Flag,
            "Indicates a novel structural variation",
        )),
        Standard::EndPosition => Some((
            Number::Count(1),
            Type::Integer,
            "End position of the variant described in this record",
        )),
        Standard::SvType => Some((Number::Count(1), Type::String, "Type of structural variant")),
        Standard::SvLengths => Some((
            Number::Unknown,
            Type::Integer,
            "Difference in length between REF and ALT alleles",
        )),
        Standard::PositionConfidenceIntervals => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around POS for imprecise variants",
        )),
        Standard::EndConfidenceIntervals => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around END for imprecise variants",
        )),
        Standard::MicrohomologyLengths => Some((
            Number::Unknown,
            Type::Integer,
            "Length of base pair identical micro-homology at event breakpoints",
        )),
        Standard::MicrohomologySequences => Some((
            Number::Unknown,
            Type::String,
            "Sequence of base pair identical micro-homology at event breakpoints",
        )),
        Standard::BreakpointIds => Some((
            Number::Unknown,
            Type::String,
            "ID of the assembled alternate allele in the assembly file",
        )),
        Standard::MobileElementInfo => Some((
            Number::Count(4),
            Type::String,
            "Mobile element info of the form NAME,START,END,POLARITY",
        )),
        Standard::MobileElementTransductionInfo => Some((
            Number::Count(4),
            Type::String,
            "Mobile element transduction info of the form CHR,START,END,POLARITY",
        )),
        Standard::DbvId => Some((
            Number::Count(1),
            Type::String,
            "ID of this element in Database of Genomic Variation",
        )),
        Standard::DbVarId => Some((
            Number::Count(1),
            Type::String,
            "ID of this element in DBVAR",
        )),
        Standard::DbRipId => Some((
            Number::Count(1),
            Type::String,
            "ID of this element in DBRIP",
        )),
        Standard::MateBreakendIds => Some((Number::Unknown, Type::String, "ID of mate breakends")),
        Standard::PartnerBreakendId => {
            Some((Number::Count(1), Type::String, "ID of partner breakend"))
        }
        Standard::BreakendEventId => Some((
            Number::Count(1),
            Type::String,
            "ID of event associated to breakend",
        )),
        Standard::BreakendConfidenceIntervals => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around the inserted material between breakends",
        )),
        Standard::AdjacentReadDepths => {
            Some((Number::Unknown, Type::Integer, "Read Depth of adjacency"))
        }
        Standard::BreakendCopyNumber => Some((
            Number::Count(1),
            Type::Integer,
            "Copy number of segment containing breakend",
        )),
        Standard::AdjacentCopyNumber => {
            Some((Number::Unknown, Type::Integer, "Copy number of adjacency"))
        }
        Standard::CopyNumberConfidenceIntervals => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around copy number for the segment",
        )),
        Standard::AdjacentCopyNumberConfidenceIntervals => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around copy number for the adjacency",
        )),

        _ => None,
    }
}
