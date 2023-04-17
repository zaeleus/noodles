use crate::{
    header::{record::value::map::format::Type, Number},
    record::genotypes::keys::key::Standard,
};

pub(super) fn definition(key: Standard) -> Option<(Number, Type, &'static str)> {
    match key {
        Standard::ReadDepths => Some((Number::R, Type::Integer, "Read depth for each allele")),
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
        Standard::ReadDepth => Some((Number::Count(1), Type::Integer, "Read depth")),
        Standard::ExpectedAlternateAlleleCounts => {
            Some((Number::A, Type::Integer, "Expected alternate allele counts"))
        }
        Standard::Filter => Some((
            Number::Count(1),
            Type::String,
            r#"Filter indicating if this genotype was "called""#,
        )),
        Standard::GenotypeLikelihoods => Some((Number::G, Type::Float, "Genotype likelihoods")),
        Standard::GenotypePosteriorProbabilities => {
            Some((Number::G, Type::Float, "Genotype posterior probabilities"))
        }
        Standard::ConditionalGenotypeQuality => Some((
            Number::Count(1),
            Type::Integer,
            "Conditional genotype quality",
        )),
        Standard::Genotype => Some((Number::Count(1), Type::String, "Genotype")),
        Standard::HaplotypeQuality => Some((Number::Count(2), Type::Integer, "Haplotype quality")),
        Standard::MappingQuality => Some((Number::Count(1), Type::Integer, "RMS mapping quality")),
        Standard::RoundedGenotypeLikelihoods => Some((
            Number::G,
            Type::Integer,
            "Phred-scaled genotype likelihoods rounded to the closest integer",
        )),
        Standard::RoundedGenotypePosteriorProbabilities => Some((
            Number::G,
            Type::Integer,
            "Phred-scaled genotype posterior probabilities rounded to the closest integer",
        )),
        Standard::PhasingQuality => Some((Number::Count(1), Type::Integer, "Phasing quality")),
        Standard::PhaseSet => Some((Number::Count(1), Type::Integer, "Phase set")),

        Standard::GenotypeCopyNumber => Some((
            Number::Count(1),
            Type::Integer,
            "Copy number genotype for imprecise events",
        )),
        Standard::GenotypeCopyNumberQuality => Some((
            Number::Count(1),
            Type::Float,
            "Copy number genotype quality for imprecise events",
        )),
        Standard::GenotypeCopyNumberLikelihoods => Some((
            Number::G,
            Type::Float,
            "Copy number genotype likelihood for imprecise events",
        )),
        Standard::GenotypeCopyNumberPosteriorProbabilities => Some((
            Number::G,
            Type::Float,
            "Copy number posterior probabilities",
        )),
        Standard::NovelVariantQualityScore => Some((
            Number::Count(1),
            Type::Integer,
            "Phred style probability score that the variant is novel",
        )),
        Standard::HaplotypeId => Some((
            Number::Count(1),
            Type::Integer,
            "Unique haplotype identifier",
        )),
        Standard::AncestralHaplotypeId => Some((
            Number::Count(1),
            Type::Integer,
            "Unique identifier of ancestral haplotype",
        )),

        _ => None,
    }
}
