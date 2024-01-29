use crate::{
    header::{record::value::map::format::Type, Number},
    record::genotypes::keys::key,
};

pub(super) fn definition(key: &str) -> Option<(Number, Type, &'static str)> {
    match key {
        key::READ_DEPTHS => Some((Number::R, Type::Integer, "Read depth for each allele")),
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
        key::READ_DEPTH => Some((Number::Count(1), Type::Integer, "Read depth")),
        key::EXPECTED_ALTERNATE_ALLELE_COUNTS => {
            Some((Number::A, Type::Integer, "Expected alternate allele counts"))
        }
        key::FILTER => Some((
            Number::Count(1),
            Type::String,
            r#"Filter indicating if this genotype was "called""#,
        )),
        key::GENOTYPE_LIKELIHOODS => Some((Number::G, Type::Float, "Genotype likelihoods")),
        key::GENOTYPE_POSTERIOR_PROBABILITIES => {
            Some((Number::G, Type::Float, "Genotype posterior probabilities"))
        }
        key::CONDITIONAL_GENOTYPE_QUALITY => Some((
            Number::Count(1),
            Type::Integer,
            "Conditional genotype quality",
        )),
        key::GENOTYPE => Some((Number::Count(1), Type::String, "Genotype")),
        key::HAPLOTYPE_QUALITY => Some((Number::Count(2), Type::Integer, "Haplotype quality")),
        key::MAPPING_QUALITY => Some((Number::Count(1), Type::Integer, "RMS mapping quality")),
        key::ROUNDED_GENOTYPE_LIKELIHOODS => Some((
            Number::G,
            Type::Integer,
            "Phred-scaled genotype likelihoods rounded to the closest integer",
        )),
        key::ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES => Some((
            Number::G,
            Type::Integer,
            "Phred-scaled genotype posterior probabilities rounded to the closest integer",
        )),
        key::PHASING_QUALITY => Some((Number::Count(1), Type::Integer, "Phasing quality")),
        key::PHASE_SET => Some((Number::Count(1), Type::Integer, "Phase set")),
        // FIXME: Number::P
        key::PHASE_SET_LIST => Some((Number::Unknown, Type::String, "Phase set list")),
        key::PHASE_SET_LIST_ORDINALS => Some((
            Number::Unknown, // FIXME: Number::P
            Type::Integer,
            "Phase set list ordinal",
        )),
        key::PHASE_SET_LIST_QUALITIES => Some((
            Number::Unknown, // FIX: Number::P
            Type::Integer,
            "Phase set list quality",
        )),

        key::GENOTYPE_COPY_NUMBER => Some((Number::Count(1), Type::Float, "Copy number")),
        key::COPY_NUMBER_CONFIDENCE_INTERVAL => Some((
            Number::Count(2),
            Type::Float,
            "Confidence interval around copy number",
        )),
        key::GENOTYPE_COPY_NUMBER_QUALITY => Some((
            Number::Count(1),
            Type::Float,
            "Copy number genotype quality",
        )),
        key::GENOTYPE_COPY_NUMBER_LIKELIHOODS => {
            Some((Number::G, Type::Float, "Copy number genotype likelihood"))
        }
        key::GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES => Some((
            Number::G,
            Type::Float,
            "Copy number posterior probabilities",
        )),
        key::NOVEL_VARIANT_QUALITY_SCORE => Some((
            Number::Count(1),
            Type::Integer,
            "Phred style probability score that the variant is novel",
        )),
        key::HAPLOTYPE_ID => Some((
            Number::Count(1),
            Type::Integer,
            "Unique haplotype identifier",
        )),
        key::ANCESTRAL_HAPLOTYPE_ID => Some((
            Number::Count(1),
            Type::Integer,
            "Unique identifier of ancestral haplotype",
        )),

        _ => None,
    }
}
