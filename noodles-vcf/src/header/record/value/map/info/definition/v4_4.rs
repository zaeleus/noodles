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
            "End position of the longest variant described in this record",
        )),
        Standard::SvType => Some((Number::Count(1), Type::String, "Type of structural variant")),
        Standard::SvLengths => Some((Number::A, Type::Integer, "Length of structural variant")),
        Standard::PositionConfidenceIntervals => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around POS for symbolic structural variants",
        )),
        Standard::EndConfidenceIntervals => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around END for symbolic structural variants",
        )),
        Standard::MicrohomologyLengths => Some((
            Number::A,
            Type::Integer,
            "Length of base pair identical micro-homology at breakpoints",
        )),
        Standard::MicrohomologySequences => Some((
            Number::A,
            Type::String,
            "Sequence of base pair identical micro-homology at breakpoints",
        )),
        Standard::BreakpointIds => Some((
            Number::A,
            Type::String,
            "ID of the assembled alternate allele in the assembly file",
        )),
        Standard::MobileElementInfo => Some((
            Number::Unknown,
            Type::String,
            "Mobile element info of the form NAME,START,END,POLARITY",
        )),
        Standard::MobileElementTransductionInfo => Some((
            Number::Unknown,
            Type::String,
            "Mobile element transduction info of the form CHR,START,END,POLARITY",
        )),
        Standard::DbvId => Some((
            Number::A,
            Type::String,
            "ID of this element in Database of Genomic Variation",
        )),
        Standard::DbVarId => Some((Number::A, Type::String, "ID of this element in DBVAR")),
        Standard::DbRipId => Some((Number::A, Type::String, "ID of this element in DBRIP")),
        Standard::MateBreakendIds => Some((Number::A, Type::String, "ID of mate breakend")),
        Standard::PartnerBreakendId => Some((Number::A, Type::String, "ID of partner breakend")),
        Standard::BreakendEventId => Some((Number::A, Type::String, "ID of associated event")),
        Standard::EventType => Some((Number::A, Type::String, "Type of associated event")),
        Standard::BreakendConfidenceIntervals => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval for the SVLEN field",
        )),
        Standard::BreakendCopyNumber => {
            Some((Number::A, Type::Float, "Copy number of CNV/breakpoint"))
        }
        Standard::CopyNumberConfidenceIntervals => Some((
            Number::Unknown,
            Type::Float,
            "Confidence interval around copy number",
        )),
        Standard::SvClaim => Some((
            Number::A,
            Type::String,
            "Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively",
        )),
        Standard::TotalRepeatSequenceCounts => Some((
            Number::A,
            Type::Integer,
            "Total number of repeat sequences in this allele",
        )),
        Standard::RepeatUnitSequences => Some((
            Number::Unknown,
            Type::String,
            "Repeat unit sequence of the corresponding repeat sequence",
        )),
        Standard::RepeatUnitLengths => Some((
            Number::Unknown,
            Type::Integer,
            "Repeat unit length of the corresponding repeat sequence",
        )),
        Standard::RepeatUnitCounts => Some((
            Number::Unknown,
            Type::Float,
            "Repeat unit count of corresponding repeat sequence",
        )),
        Standard::TotalRepeatSequenceBaseCounts => Some((
            Number::Unknown,
            Type::Integer,
            "Total number of bases in the corresponding repeat sequence",
        )),
        Standard::RepeatUnitCountConfidenceIntervals => Some((
            Number::Unknown,
            Type::Float,
            "Confidence interval around RUC",
        )),
        Standard::TotalRepeatSequenceBaseCountConfidenceIntervals => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around RB",
        )),
        Standard::RepeatUnitBaseCounts => Some((
            Number::Unknown,
            Type::Integer,
            "Number of bases in each individual repeat unit",
        )),

        _ => None,
    }
}
