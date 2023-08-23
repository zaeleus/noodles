use std::str::FromStr;

use super::ParseError;

/// A reserved VCF header format key.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Standard {
    /// Read depth for each allele (`AD`).
    ReadDepths,
    /// Read depth for each allele on the forward strand (`ADF`).
    ForwardStrandReadDepths,
    /// Read depth for each allele on the reverse strand (`ADR`).
    ReverseStrandReadDepths,
    /// Read depth (`DP`).
    ReadDepth,
    /// Expected alternate allele counts (`EC`).
    ExpectedAlternateAlleleCounts,
    /// Filter indicating if this genotype was "called" (`FT`).
    Filter,
    /// Genotype likelihoods (`GL`).
    GenotypeLikelihoods,
    /// Genotype posterior probabilities (`GP`).
    GenotypePosteriorProbabilities,
    /// Conditional genotype quality (`GQ`).
    ConditionalGenotypeQuality,
    /// Genotype (`GT`).
    Genotype,
    /// Haplotype quality (`HQ`).
    HaplotypeQuality,
    /// RMS mapping quality (`MQ`).
    MappingQuality,
    /// Phred-scaled genotype likelihoods rounded to the closest integer (`PL`).
    RoundedGenotypeLikelihoods,
    /// Phred-scaled genotype posterior probabilities rounded to the closest integer (`PP`).
    RoundedGenotypePosteriorProbabilities,
    /// Phasing quality (`PQ`).
    PhasingQuality,
    /// Phase set (`PS`).
    PhaseSet,
    /// Phase set list (`PSL`).
    ///
    /// Added in VCF 4.4.
    PhaseSetList,
    /// Phase set list ordinal (`PSO`).
    ///
    /// Added in VCF 4.4.
    PhaseSetListOrdinals,
    /// Phase set list quality (`PSQ`).
    ///
    /// Added in VCF 4.4.
    PhaseSetListQualities,

    /// Copy number genotype for imprecise events (`CN`).
    GenotypeCopyNumber,
    /// Confidence interval around copy number (`CICN`).
    ///
    /// Added in VCF 4.4.
    CopyNumberConfidenceInterval,
    /// Copy number genotype quality for imprecise events (`CNQ`).
    GenotypeCopyNumberQuality,
    /// Copy number genotype likelihood for imprecise events (`CNL`).
    GenotypeCopyNumberLikelihoods,
    /// Copy number posterior probabilities (`CNP`).
    GenotypeCopyNumberPosteriorProbabilities,
    /// Phred style probability score that the variant is novel (`NQ`).
    NovelVariantQualityScore,
    /// Unique haplotype identifier (`HAP`).
    HaplotypeId,
    /// Unique identifier of ancestral haplotype (`AHAP`).
    AncestralHaplotypeId,
}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::ReadDepths => "AD",
            Self::ForwardStrandReadDepths => "ADF",
            Self::ReverseStrandReadDepths => "ADR",
            Self::ReadDepth => "DP",
            Self::ExpectedAlternateAlleleCounts => "EC",
            Self::Filter => "FT",
            Self::GenotypeLikelihoods => "GL",
            Self::GenotypePosteriorProbabilities => "GP",
            Self::ConditionalGenotypeQuality => "GQ",
            Self::Genotype => "GT",
            Self::HaplotypeQuality => "HQ",
            Self::MappingQuality => "MQ",
            Self::RoundedGenotypeLikelihoods => "PL",
            Self::RoundedGenotypePosteriorProbabilities => "PP",
            Self::PhasingQuality => "PQ",
            Self::PhaseSet => "PS",
            Self::PhaseSetList => "PSL",
            Self::PhaseSetListOrdinals => "PSO",
            Self::PhaseSetListQualities => "PSQ",

            Self::GenotypeCopyNumber => "CN",
            Self::CopyNumberConfidenceInterval => "CICN",
            Self::GenotypeCopyNumberQuality => "CNQ",
            Self::GenotypeCopyNumberLikelihoods => "CNL",
            Self::GenotypeCopyNumberPosteriorProbabilities => "CNP",
            Self::NovelVariantQualityScore => "NQ",
            Self::HaplotypeId => "HAP",
            Self::AncestralHaplotypeId => "AHAP",
        }
    }
}

impl FromStr for Standard {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "AD" => Ok(Self::ReadDepths),
            "ADF" => Ok(Self::ForwardStrandReadDepths),
            "ADR" => Ok(Self::ReverseStrandReadDepths),
            "DP" => Ok(Self::ReadDepth),
            "EC" => Ok(Self::ExpectedAlternateAlleleCounts),
            "FT" => Ok(Self::Filter),
            "GL" => Ok(Self::GenotypeLikelihoods),
            "GP" => Ok(Self::GenotypePosteriorProbabilities),
            "GQ" => Ok(Self::ConditionalGenotypeQuality),
            "GT" => Ok(Self::Genotype),
            "HQ" => Ok(Self::HaplotypeQuality),
            "MQ" => Ok(Self::MappingQuality),
            "PL" => Ok(Self::RoundedGenotypeLikelihoods),
            "PP" => Ok(Self::RoundedGenotypePosteriorProbabilities),
            "PQ" => Ok(Self::PhasingQuality),
            "PS" => Ok(Self::PhaseSet),
            "PSL" => Ok(Self::PhaseSetList),
            "PSO" => Ok(Self::PhaseSetListOrdinals),
            "PSQ" => Ok(Self::PhaseSetListQualities),

            "CN" => Ok(Self::GenotypeCopyNumber),
            "CICN" => Ok(Self::CopyNumberConfidenceInterval),
            "CNQ" => Ok(Self::GenotypeCopyNumberQuality),
            "CNL" => Ok(Self::GenotypeCopyNumberLikelihoods),
            "CNP" => Ok(Self::GenotypeCopyNumberPosteriorProbabilities),
            "NQ" => Ok(Self::NovelVariantQualityScore),
            "HAP" => Ok(Self::HaplotypeId),
            "AHAP" => Ok(Self::AncestralHaplotypeId),

            _ => Err(ParseError::Invalid),
        }
    }
}
