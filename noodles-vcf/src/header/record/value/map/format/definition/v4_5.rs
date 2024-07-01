use crate::{
    header::record::value::map::format::{Number, Type},
    variant::record::samples::keys::key,
};

pub(super) fn definition(key: &str) -> Option<(Number, Type, &'static str)> {
    match key {
        key::READ_DEPTHS => Some((
            Number::ReferenceAlternateBases,
            Type::Integer,
            "Read depth for each allele",
        )),
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
        key::READ_DEPTH => Some((Number::Count(1), Type::Integer, "Read depth")),
        key::EXPECTED_ALTERNATE_ALLELE_COUNTS => Some((
            Number::AlternateBases,
            Type::Integer,
            "Expected alternate allele counts",
        )),
        key::LENGTH => Some((
            Number::Count(1),
            Type::Integer,
            "Length of <*> reference block",
        )),
        key::FILTER => Some((
            Number::Count(1),
            Type::String,
            r#"Filter indicating if this genotype was "called""#,
        )),
        key::GENOTYPE_LIKELIHOODS => Some((Number::Samples, Type::Float, "Genotype likelihoods")),
        key::GENOTYPE_POSTERIOR_PROBABILITIES => Some((
            Number::Samples,
            Type::Float,
            "Genotype posterior probabilities",
        )),
        key::CONDITIONAL_GENOTYPE_QUALITY => Some((
            Number::Count(1),
            Type::Integer,
            "Conditional genotype quality",
        )),
        key::GENOTYPE => Some((Number::Count(1), Type::String, "Genotype")),
        key::HAPLOTYPE_QUALITY => Some((Number::Count(2), Type::Integer, "Haplotype quality")),
        key::RESERVED_LA => Some((Number::Unknown, Type::Integer, "Reserved")),
        key::LOCAL_ALTERNATIVE_ALLELE => Some((
            Number::Unknown,
            Type::Integer,
            "1-based indices into ALT, indicating which alleles are relevant (local) for the current sample",
        )),
        key::LOCAL_READ_DEPTHS => Some((
            Number::LocalReferenceAlternateBases,
            Type::Integer,
            "Local-allele representation of AD",
        )),
        key::LOCAL_FORWARD_STRAND_READ_DEPTHS => Some((
            Number::LocalReferenceAlternateBases,
            Type::Integer,
            "Local-allele representation of ADF",
        )),
        key::LOCAL_REVERSE_STRAND_READ_DEPTHS => Some((
            Number::LocalReferenceAlternateBases,
            Type::Integer,
            "Local-allele representation of ADR",
        )),
        key::LOCAL_EXPECTED_ALTERNATE_ALLELE_COUNTS => Some((
            Number::LocalAlternateBases,
            Type::Integer,
            "Local-allele representation of EC",
        )),
        key::LOCAL_GENOTYPE_LIKELIHOODS => Some((
            Number::LocalSamples,
            Type::Integer,
            "Local-allele representation of GL",
        )),
        key::LOCAL_GENOTYPE_POSTERIOR_PROBABILITIES => Some((
            Number::LocalSamples,
            Type::Integer,
            "Local-allele representation of GP",
        )),
        key::LOCAL_ROUNDED_GENOTYPE_LIKELIHOODS => Some((
            Number::LocalSamples,
            Type::Integer,
            "Local-allele representation of PL",
        )),
        key::LOCAL_ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES => Some((
            Number::LocalSamples,
            Type::Integer,
            "Local-allele representation of PP",
        )),
        key::MAPPING_QUALITY => Some((Number::Count(1), Type::Integer, "RMS mapping quality")),
        key::ROUNDED_GENOTYPE_LIKELIHOODS => Some((
            Number::Samples,
            Type::Integer,
            "Phred-scaled genotype likelihoods rounded to the closest integer",
        )),
        key::ROUNDED_GENOTYPE_POSTERIOR_PROBABILITIES => Some((
            Number::Samples,
            Type::Integer,
            "Phred-scaled genotype posterior probabilities rounded to the closest integer",
        )),
        key::PHASING_QUALITY => Some((Number::Count(1), Type::Integer, "Phasing quality")),
        key::PHASE_SET => Some((Number::Count(1), Type::Integer, "Phase set")),
        key::PHASE_SET_LIST => Some((Number::Ploidy, Type::String, "Phase set list")),
        key::PHASE_SET_LIST_ORDINALS => Some((
            Number::Ploidy,
            Type::Integer,
            "Phase set list ordinal",
        )),
        key::PHASE_SET_LIST_QUALITIES => Some((
            Number::Ploidy,
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
        key::GENOTYPE_COPY_NUMBER_LIKELIHOODS => Some((
            Number::Samples,
            Type::Float,
            "Copy number genotype likelihood",
        )),
        key::GENOTYPE_COPY_NUMBER_POSTERIOR_PROBABILITIES => Some((
            Number::Samples,
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
