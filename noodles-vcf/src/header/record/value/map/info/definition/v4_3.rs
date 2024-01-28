use crate::{
    header::{record::value::map::info::Type, Number},
    record::info::field::key,
};

pub(super) fn definition(key: &str) -> Option<(Number, Type, &'static str)> {
    match key {
        key::ANCESTRAL_ALLELE => Some((Number::Count(1), Type::String, "Ancestral allele")),
        key::ALLELE_COUNT => Some((
            Number::A,
            Type::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        )),
        key::TOTAL_READ_DEPTHS => {
            Some((Number::R, Type::Integer, "Total read depth for each allele"))
        }
        key::FORWARD_STRAND_READ_DEPTHS => Some((
            Number::R,
            Type::Integer,
            "Read depth for each allele on the forward strand",
        )),
        key::REVERSE_STRAND_READ_DEPTHS => Some((
            Number::R,
            Type::Integer,
            "Read depth for each allele on the reverse strand",
        )),
        key::ALLELE_FREQUENCIES => Some((
            Number::A,
            Type::Float,
            "Allele frequency for each ALT allele in the same order as listed",
        )),
        key::TOTAL_ALLELE_COUNT => Some((
            Number::Count(1),
            Type::Integer,
            "Total number of alleles in called genotypes",
        )),
        key::BASE_QUALITY => Some((Number::Count(1), Type::Float, "RMS base quality")),
        key::CIGAR => Some((
            Number::A,
            Type::String,
            "Cigar string describing how to align an alternate allele to the reference allele",
        )),
        key::IS_IN_DB_SNP => Some((Number::Count(0), Type::Flag, "dbSNP membership")),
        key::TOTAL_DEPTH => Some((
            Number::Count(1),
            Type::Integer,
            "Combined depth across samples",
        )),
        key::IS_IN_HAP_MAP_2 => Some((Number::Count(0), Type::Flag, "HapMap2 membership")),
        key::IS_IN_HAP_MAP_3 => Some((Number::Count(0), Type::Flag, "HapMap3 membership")),
        key::MAPPING_QUALITY => Some((Number::Count(1), Type::Float, "RMS mapping quality")),
        key::ZERO_MAPPING_QUALITY_COUNT => {
            Some((Number::Count(1), Type::Integer, "Number of MAPQ == 0 reads"))
        }
        key::SAMPLES_WITH_DATA_COUNT => Some((
            Number::Count(1),
            Type::Integer,
            "Number of samples with data",
        )),
        key::STRAND_BIAS => Some((Number::Count(4), Type::Integer, "Strand bias")),
        key::IS_SOMATIC_MUTATION => Some((Number::Count(0), Type::Flag, "Somatic mutation")),
        key::IS_VALIDATED => Some((
            Number::Count(0),
            Type::Flag,
            "Validated by follow-up experiment",
        )),
        key::IS_IN_1000_GENOMES => Some((Number::Count(0), Type::Flag, "1000 Genomes membership")),

        key::IS_IMPRECISE => Some((
            Number::Count(0),
            Type::Flag,
            "Imprecise structural variation",
        )),
        key::IS_NOVEL => Some((
            Number::Count(0),
            Type::Flag,
            "Indicates a novel structural variation",
        )),
        key::END_POSITION => Some((
            Number::Count(1),
            Type::Integer,
            "End position of the variant described in this record",
        )),
        key::SV_TYPE => Some((Number::Count(1), Type::String, "Type of structural variant")),
        key::SV_LENGTHS => Some((
            Number::Unknown,
            Type::Integer,
            "Difference in length between REF and ALT alleles",
        )),
        key::POSITION_CONFIDENCE_INTERVALS => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around POS for imprecise variants",
        )),
        key::END_CONFIDENCE_INTERVALS => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around END for imprecise variants",
        )),
        key::MICROHOMOLOGY_LENGTHS => Some((
            Number::Unknown,
            Type::Integer,
            "Length of base pair identical micro-homology at event breakpoints",
        )),
        key::MICROHOMOLOGY_SEQUENCES => Some((
            Number::Unknown,
            Type::String,
            "Sequence of base pair identical micro-homology at event breakpoints",
        )),
        key::BREAKPOINT_IDS => Some((
            Number::Unknown,
            Type::String,
            "ID of the assembled alternate allele in the assembly file",
        )),
        key::MOBILE_ELEMENT_INFO => Some((
            Number::Count(4),
            Type::String,
            "Mobile element info of the form NAME,START,END,POLARITY",
        )),
        key::MOBILE_ELEMENT_TRANSDUCTION_INFO => Some((
            Number::Count(4),
            Type::String,
            "Mobile element transduction info of the form CHR,START,END,POLARITY",
        )),
        key::DBV_ID => Some((
            Number::Count(1),
            Type::String,
            "ID of this element in Database of Genomic Variation",
        )),
        key::DB_VAR_ID => Some((
            Number::Count(1),
            Type::String,
            "ID of this element in DBVAR",
        )),
        key::DB_RIP_ID => Some((
            Number::Count(1),
            Type::String,
            "ID of this element in DBRIP",
        )),
        key::MATE_BREAKEND_IDS => Some((Number::Unknown, Type::String, "ID of mate breakends")),
        key::PARTNER_BREAKEND_ID => {
            Some((Number::Count(1), Type::String, "ID of partner breakend"))
        }
        key::BREAKEND_EVENT_ID => Some((
            Number::Count(1),
            Type::String,
            "ID of event associated to breakend",
        )),
        key::BREAKEND_CONFIDENCE_INTERVALS => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around the inserted material between breakends",
        )),
        key::ADJACENT_READ_DEPTHS => {
            Some((Number::Unknown, Type::Integer, "Read Depth of adjacency"))
        }
        key::BREAKEND_COPY_NUMBER => Some((
            Number::Count(1),
            Type::Integer,
            "Copy number of segment containing breakend",
        )),
        key::ADJACENT_COPY_NUMBER => {
            Some((Number::Unknown, Type::Integer, "Copy number of adjacency"))
        }
        key::COPY_NUMBER_CONFIDENCE_INTERVALS => Some((
            Number::Count(2),
            Type::Integer,
            "Confidence interval around copy number for the segment",
        )),
        key::ADJACENT_COPY_NUMBER_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around copy number for the adjacency",
        )),

        _ => None,
    }
}
