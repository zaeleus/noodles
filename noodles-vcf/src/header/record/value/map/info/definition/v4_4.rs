use crate::{
    header::record::value::map::info::{Number, Type},
    variant::record::info::field::key,
};

pub(super) fn definition(key: &str) -> Option<(Number, Type, &'static str)> {
    match key {
        key::ANCESTRAL_ALLELE => Some((Number::Count(1), Type::String, "Ancestral allele")),
        key::ALLELE_COUNT => Some((
            Number::AlternateBases,
            Type::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        )),
        key::TOTAL_READ_DEPTHS => {
            Some((Number::ReferenceAlternateBases, Type::Integer, "Total read depth for each allele"))
        }
        key::FORWARD_STRAND_READ_DEPTHS => Some((
            Number::ReferenceAlternateBases,
            Type::Integer,
            "Read depth for each allele on the forward strand",
        )),
        key::REVERSE_STRAND_READ_DEPTHS => Some((
            Number::ReferenceAlternateBases,
            Type::Integer,
            "Read depth for each allele on the reverse strand",
        )),
        key::ALLELE_FREQUENCIES => Some((
            Number::AlternateBases,
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
            Number::AlternateBases,
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
        key::IS_IN_1000_GENOMES => {
            Some((Number::Count(0), Type::Flag, "1000 Genomes membership"))
        }

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
            "End position of the longest variant described in this record",
        )),
        key::SV_TYPE => Some((Number::Count(1), Type::String, "Type of structural variant")),
        key::SV_LENGTHS => Some((Number::AlternateBases, Type::Integer, "Length of structural variant")),
        key::POSITION_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around POS for symbolic structural variants",
        )),
        key::END_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around END for symbolic structural variants",
        )),
        key::MICROHOMOLOGY_LENGTHS => Some((
            Number::AlternateBases,
            Type::Integer,
            "Length of base pair identical micro-homology at breakpoints",
        )),
        key::MICROHOMOLOGY_SEQUENCES => Some((
            Number::AlternateBases,
            Type::String,
            "Sequence of base pair identical micro-homology at breakpoints",
        )),
        key::BREAKPOINT_IDS => Some((
            Number::AlternateBases,
            Type::String,
            "ID of the assembled alternate allele in the assembly file",
        )),
        key::MOBILE_ELEMENT_INFO => Some((
            Number::Unknown,
            Type::String,
            "Mobile element info of the form NAME,START,END,POLARITY",
        )),
        key::MOBILE_ELEMENT_TRANSDUCTION_INFO => Some((
            Number::Unknown,
            Type::String,
            "Mobile element transduction info of the form CHR,START,END,POLARITY",
        )),
        key::DBV_ID => Some((
            Number::AlternateBases,
            Type::String,
            "ID of this element in Database of Genomic Variation",
        )),
        key::DB_VAR_ID => Some((Number::AlternateBases, Type::String, "ID of this element in DBVAR")),
        key::DB_RIP_ID => Some((Number::AlternateBases, Type::String, "ID of this element in DBRIP")),
        key::MATE_BREAKEND_IDS => Some((Number::AlternateBases, Type::String, "ID of mate breakend")),
        key::PARTNER_BREAKEND_ID => Some((Number::AlternateBases, Type::String, "ID of partner breakend")),
        key::BREAKEND_EVENT_ID => Some((Number::AlternateBases, Type::String, "ID of associated event")),
        key::EVENT_TYPE => Some((Number::AlternateBases, Type::String, "Type of associated event")),
        key::BREAKEND_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval for the SVLEN field",
        )),
        key::BREAKEND_COPY_NUMBER => {
            Some((Number::AlternateBases, Type::Float, "Copy number of CNV/breakpoint"))
        }
        key::COPY_NUMBER_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Float,
            "Confidence interval around copy number",
        )),
        key::SV_CLAIM => Some((
            Number::AlternateBases,
            Type::String,
            "Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively",
        )),
        key::TOTAL_REPEAT_SEQUENCE_COUNTS => Some((
            Number::AlternateBases,
            Type::Integer,
            "Total number of repeat sequences in this allele",
        )),
        key::REPEAT_UNIT_SEQUENCES => Some((
            Number::Unknown,
            Type::String,
            "Repeat unit sequence of the corresponding repeat sequence",
        )),
        key::REPEAT_UNIT_LENGTHS => Some((
            Number::Unknown,
            Type::Integer,
            "Repeat unit length of the corresponding repeat sequence",
        )),
        key::REPEAT_UNIT_COUNTS => Some((
            Number::Unknown,
            Type::Float,
            "Repeat unit count of corresponding repeat sequence",
        )),
        key::TOTAL_REPEAT_SEQUENCE_BASE_COUNTS => Some((
            Number::Unknown,
            Type::Integer,
            "Total number of bases in the corresponding repeat sequence",
        )),
        key::REPEAT_UNIT_COUNT_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Float,
            "Confidence interval around RUC",
        )),
        key::TOTAL_REPEAT_SEQUENCE_BASE_COUNT_CONFIDENCE_INTERVALS => Some((
            Number::Unknown,
            Type::Integer,
            "Confidence interval around RB",
        )),
        key::REPEAT_UNIT_BASE_COUNTS => Some((
            Number::Unknown,
            Type::Integer,
            "Number of bases in each individual repeat unit",
        )),

        _ => None,
    }
}
