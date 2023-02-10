//! VCF header format key.

use std::{error, fmt, str::FromStr};

use crate::header::{record::value::map::format::Type, Number};

/// A VCF header format key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Key {
    // § 1.6.2 Genotype fields (2021-01-13)
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

    // § 4 FORMAT keys used for structural variants (2021-01-13)
    /// Copy number genotype for imprecise events (`CN`).
    GenotypeCopyNumber,
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

    /// Any other non-reserved key.
    Other(String),
}

pub(crate) fn number(key: &Key) -> Option<Number> {
    match key {
        Key::ReadDepths => Some(Number::R),
        Key::ForwardStrandReadDepths => Some(Number::R),
        Key::ReverseStrandReadDepths => Some(Number::R),
        Key::ReadDepth => Some(Number::Count(1)),
        Key::ExpectedAlternateAlleleCounts => Some(Number::A),
        Key::Filter => Some(Number::Count(1)),
        Key::GenotypeLikelihoods => Some(Number::G),
        Key::GenotypePosteriorProbabilities => Some(Number::G),
        Key::ConditionalGenotypeQuality => Some(Number::Count(1)),
        Key::Genotype => Some(Number::Count(1)),
        Key::HaplotypeQuality => Some(Number::Count(2)),
        Key::MappingQuality => Some(Number::Count(1)),
        Key::RoundedGenotypeLikelihoods => Some(Number::G),
        Key::RoundedGenotypePosteriorProbabilities => Some(Number::G),
        Key::PhasingQuality => Some(Number::Count(1)),
        Key::PhaseSet => Some(Number::Count(1)),

        Key::GenotypeCopyNumber => Some(Number::Count(1)),
        Key::GenotypeCopyNumberQuality => Some(Number::Count(1)),
        Key::GenotypeCopyNumberLikelihoods => Some(Number::G),
        Key::GenotypeCopyNumberPosteriorProbabilities => Some(Number::G),
        Key::NovelVariantQualityScore => Some(Number::Count(1)),
        Key::HaplotypeId => Some(Number::Count(1)),
        Key::AncestralHaplotypeId => Some(Number::Count(1)),

        Key::Other(_) => None,
    }
}

pub(crate) fn ty(key: &Key) -> Option<Type> {
    match key {
        Key::ReadDepths => Some(Type::Integer),
        Key::ForwardStrandReadDepths => Some(Type::Integer),
        Key::ReverseStrandReadDepths => Some(Type::Integer),
        Key::ReadDepth => Some(Type::Integer),
        Key::ExpectedAlternateAlleleCounts => Some(Type::Integer),
        Key::Filter => Some(Type::String),
        Key::GenotypeLikelihoods => Some(Type::Float),
        Key::GenotypePosteriorProbabilities => Some(Type::Float),
        Key::ConditionalGenotypeQuality => Some(Type::Integer),
        Key::Genotype => Some(Type::String),
        Key::HaplotypeQuality => Some(Type::Integer),
        Key::MappingQuality => Some(Type::Integer),
        Key::RoundedGenotypeLikelihoods => Some(Type::Integer),
        Key::RoundedGenotypePosteriorProbabilities => Some(Type::Integer),
        Key::PhasingQuality => Some(Type::Integer),
        Key::PhaseSet => Some(Type::Integer),

        Key::GenotypeCopyNumber => Some(Type::Integer),
        Key::GenotypeCopyNumberQuality => Some(Type::Float),
        Key::GenotypeCopyNumberLikelihoods => Some(Type::Float),
        Key::GenotypeCopyNumberPosteriorProbabilities => Some(Type::Float),
        Key::NovelVariantQualityScore => Some(Type::Integer),
        Key::HaplotypeId => Some(Type::Integer),
        Key::AncestralHaplotypeId => Some(Type::Integer),

        Key::Other(_) => None,
    }
}

pub(crate) fn description(key: &Key) -> Option<&str> {
    match key {
        Key::ReadDepths => Some("Read depth for each allele"),
        Key::ForwardStrandReadDepths => Some("Read depth for each allele on the forward strand"),
        Key::ReverseStrandReadDepths => Some("Read depth for each allele on the reverse strand"),
        Key::ReadDepth => Some("Read depth"),
        Key::ExpectedAlternateAlleleCounts => Some("Expected alternate allele counts"),
        Key::Filter => Some(r#"Filter indicating if this genotype was "called""#),
        Key::GenotypeLikelihoods => Some("Genotype likelihoods"),
        Key::GenotypePosteriorProbabilities => Some("Genotype posterior probabilities"),
        Key::ConditionalGenotypeQuality => Some("Conditional genotype quality"),
        Key::Genotype => Some("Genotype"),
        Key::HaplotypeQuality => Some("Haplotype quality"),
        Key::MappingQuality => Some("RMS mapping quality"),
        Key::RoundedGenotypeLikelihoods => {
            Some("Phred-scaled genotype likelihoods rounded to the closest integer")
        }
        Key::RoundedGenotypePosteriorProbabilities => {
            Some("Phred-scaled genotype posterior probabilities rounded to the closest integer")
        }
        Key::PhasingQuality => Some("Phasing quality"),
        Key::PhaseSet => Some("Phase set"),

        Key::GenotypeCopyNumber => Some("Copy number genotype for imprecise events"),
        Key::GenotypeCopyNumberQuality => Some("Copy number genotype quality for imprecise events"),
        Key::GenotypeCopyNumberLikelihoods => {
            Some("Copy number genotype likelihood for imprecise events")
        }
        Key::GenotypeCopyNumberPosteriorProbabilities => {
            Some("Copy number posterior probabilities")
        }
        Key::NovelVariantQualityScore => {
            Some("Phred style probability score that the variant is novel")
        }
        Key::HaplotypeId => Some("Unique haplotype identifier"),
        Key::AncestralHaplotypeId => Some("Unique identifier of ancestral haplotype"),

        Key::Other(_) => None,
    }
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Key::ReadDepths => "AD",
            Key::ForwardStrandReadDepths => "ADF",
            Key::ReverseStrandReadDepths => "ADR",
            Key::ReadDepth => "DP",
            Key::ExpectedAlternateAlleleCounts => "EC",
            Key::Filter => "FT",
            Key::GenotypeLikelihoods => "GL",
            Key::GenotypePosteriorProbabilities => "GP",
            Key::ConditionalGenotypeQuality => "GQ",
            Key::Genotype => "GT",
            Key::HaplotypeQuality => "HQ",
            Key::MappingQuality => "MQ",
            Key::RoundedGenotypeLikelihoods => "PL",
            Key::RoundedGenotypePosteriorProbabilities => "PP",
            Key::PhasingQuality => "PQ",
            Key::PhaseSet => "PS",

            Key::GenotypeCopyNumber => "CN",
            Key::GenotypeCopyNumberQuality => "CNQ",
            Key::GenotypeCopyNumberLikelihoods => "CNL",
            Key::GenotypeCopyNumberPosteriorProbabilities => "CNP",
            Key::NovelVariantQualityScore => "NQ",
            Key::HaplotypeId => "HAP",
            Key::AncestralHaplotypeId => "AHAP",

            Key::Other(s) => s,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF record genotype field key fails to parse.
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

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        match s {
            "AD" => Ok(Key::ReadDepths),
            "ADF" => Ok(Key::ForwardStrandReadDepths),
            "ADR" => Ok(Key::ReverseStrandReadDepths),
            "DP" => Ok(Key::ReadDepth),
            "EC" => Ok(Key::ExpectedAlternateAlleleCounts),
            "FT" => Ok(Key::Filter),
            "GL" => Ok(Key::GenotypeLikelihoods),
            "GP" => Ok(Key::GenotypePosteriorProbabilities),
            "GQ" => Ok(Key::ConditionalGenotypeQuality),
            "GT" => Ok(Key::Genotype),
            "HQ" => Ok(Key::HaplotypeQuality),
            "MQ" => Ok(Key::MappingQuality),
            "PL" => Ok(Key::RoundedGenotypeLikelihoods),
            "PP" => Ok(Key::RoundedGenotypePosteriorProbabilities),
            "PQ" => Ok(Key::PhasingQuality),
            "PS" => Ok(Key::PhaseSet),

            "CN" => Ok(Key::GenotypeCopyNumber),
            "CNQ" => Ok(Key::GenotypeCopyNumberQuality),
            "CNL" => Ok(Key::GenotypeCopyNumberLikelihoods),
            "CNP" => Ok(Key::GenotypeCopyNumberPosteriorProbabilities),
            "NQ" => Ok(Key::NovelVariantQualityScore),
            "HAP" => Ok(Key::HaplotypeId),
            "AHAP" => Ok(Key::AncestralHaplotypeId),

            _ => {
                if is_valid_name(s) {
                    Ok(Key::Other(s.into()))
                } else {
                    Err(ParseError::Invalid)
                }
            }
        }
    }
}

// § 1.6.2 Genotype fields
fn is_valid_name_char(c: char) -> bool {
    matches!(c, '0'..='9' | 'A'..='Z' | 'a'..='z' | '_' | '.')
}

fn is_valid_name(s: &str) -> bool {
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if !matches!(c, 'A'..='Z' | 'a'..='z' | '_') {
            return false;
        }
    }

    chars.all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::ReadDepths.to_string(), "AD");
        assert_eq!(Key::ForwardStrandReadDepths.to_string(), "ADF");
        assert_eq!(Key::ReverseStrandReadDepths.to_string(), "ADR");
        assert_eq!(Key::ReadDepth.to_string(), "DP");
        assert_eq!(Key::ExpectedAlternateAlleleCounts.to_string(), "EC");
        assert_eq!(Key::Filter.to_string(), "FT");
        assert_eq!(Key::GenotypeLikelihoods.to_string(), "GL");
        assert_eq!(Key::GenotypePosteriorProbabilities.to_string(), "GP");
        assert_eq!(Key::ConditionalGenotypeQuality.to_string(), "GQ");
        assert_eq!(Key::Genotype.to_string(), "GT");
        assert_eq!(Key::HaplotypeQuality.to_string(), "HQ");
        assert_eq!(Key::MappingQuality.to_string(), "MQ");
        assert_eq!(Key::RoundedGenotypeLikelihoods.to_string(), "PL");
        assert_eq!(Key::RoundedGenotypePosteriorProbabilities.to_string(), "PP");
        assert_eq!(Key::PhasingQuality.to_string(), "PQ");
        assert_eq!(Key::PhaseSet.to_string(), "PS");

        assert_eq!(Key::GenotypeCopyNumber.to_string(), "CN");
        assert_eq!(Key::GenotypeCopyNumberQuality.to_string(), "CNQ");
        assert_eq!(Key::GenotypeCopyNumberLikelihoods.to_string(), "CNL");
        assert_eq!(
            Key::GenotypeCopyNumberPosteriorProbabilities.to_string(),
            "CNP"
        );
        assert_eq!(Key::NovelVariantQualityScore.to_string(), "NQ");
        assert_eq!(Key::HaplotypeId.to_string(), "HAP");
        assert_eq!(Key::AncestralHaplotypeId.to_string(), "AHAP");

        assert_eq!(Key::Other(String::from("NDLS")).to_string(), "NDLS");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("AD".parse(), Ok(Key::ReadDepths));
        assert_eq!("ADF".parse(), Ok(Key::ForwardStrandReadDepths));
        assert_eq!("ADR".parse(), Ok(Key::ReverseStrandReadDepths));
        assert_eq!("DP".parse(), Ok(Key::ReadDepth));
        assert_eq!("EC".parse(), Ok(Key::ExpectedAlternateAlleleCounts));
        assert_eq!("FT".parse(), Ok(Key::Filter));
        assert_eq!("GL".parse(), Ok(Key::GenotypeLikelihoods));
        assert_eq!("GP".parse(), Ok(Key::GenotypePosteriorProbabilities));
        assert_eq!("GQ".parse(), Ok(Key::ConditionalGenotypeQuality));
        assert_eq!("GT".parse(), Ok(Key::Genotype));
        assert_eq!("HQ".parse(), Ok(Key::HaplotypeQuality));
        assert_eq!("MQ".parse(), Ok(Key::MappingQuality));
        assert_eq!("PL".parse(), Ok(Key::RoundedGenotypeLikelihoods));
        assert_eq!("PP".parse(), Ok(Key::RoundedGenotypePosteriorProbabilities));
        assert_eq!("PQ".parse(), Ok(Key::PhasingQuality));
        assert_eq!("PS".parse(), Ok(Key::PhaseSet));

        assert_eq!("CN".parse(), Ok(Key::GenotypeCopyNumber));
        assert_eq!("CNQ".parse(), Ok(Key::GenotypeCopyNumberQuality));
        assert_eq!("CNL".parse(), Ok(Key::GenotypeCopyNumberLikelihoods));
        assert_eq!(
            "CNP".parse(),
            Ok(Key::GenotypeCopyNumberPosteriorProbabilities)
        );
        assert_eq!("NQ".parse(), Ok(Key::NovelVariantQualityScore));
        assert_eq!("HAP".parse(), Ok(Key::HaplotypeId));
        assert_eq!("AHAP".parse(), Ok(Key::AncestralHaplotypeId));

        assert_eq!("NDLS".parse(), Ok(Key::Other(String::from("NDLS"))));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
        assert_eq!("8D".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!(".N".parse::<Key>(), Err(ParseError::Invalid));
        assert_eq!("A!".parse::<Key>(), Err(ParseError::Invalid));
    }

    #[test]
    fn test_number() {
        assert_eq!(number(&Key::ReadDepths), Some(Number::R));
        assert_eq!(number(&Key::ForwardStrandReadDepths), Some(Number::R));
        assert_eq!(number(&Key::ReverseStrandReadDepths), Some(Number::R));
        assert_eq!(number(&Key::ReadDepth), Some(Number::Count(1)));
        assert_eq!(number(&Key::ExpectedAlternateAlleleCounts), Some(Number::A));
        assert_eq!(number(&Key::Filter), Some(Number::Count(1)));
        assert_eq!(number(&Key::GenotypeLikelihoods), Some(Number::G));
        assert_eq!(
            number(&Key::GenotypePosteriorProbabilities),
            Some(Number::G)
        );
        assert_eq!(
            number(&Key::ConditionalGenotypeQuality),
            Some(Number::Count(1))
        );
        assert_eq!(number(&Key::Genotype), Some(Number::Count(1)));
        assert_eq!(number(&Key::HaplotypeQuality), Some(Number::Count(2)));
        assert_eq!(number(&Key::MappingQuality), Some(Number::Count(1)));
        assert_eq!(number(&Key::RoundedGenotypeLikelihoods), Some(Number::G));
        assert_eq!(
            number(&Key::RoundedGenotypePosteriorProbabilities),
            Some(Number::G)
        );
        assert_eq!(number(&Key::PhasingQuality), Some(Number::Count(1)));
        assert_eq!(number(&Key::PhaseSet), Some(Number::Count(1)));

        assert_eq!(number(&Key::GenotypeCopyNumber), Some(Number::Count(1)));
        assert_eq!(
            number(&Key::GenotypeCopyNumberQuality),
            Some(Number::Count(1))
        );
        assert_eq!(number(&Key::GenotypeCopyNumberLikelihoods), Some(Number::G));
        assert_eq!(
            number(&Key::GenotypeCopyNumberPosteriorProbabilities),
            Some(Number::G)
        );
        assert_eq!(
            number(&Key::NovelVariantQualityScore),
            Some(Number::Count(1))
        );
        assert_eq!(number(&Key::HaplotypeId), Some(Number::Count(1)));
        assert_eq!(number(&Key::AncestralHaplotypeId), Some(Number::Count(1)));

        assert!(number(&Key::Other(String::from("NDLS"))).is_none());
    }

    #[test]
    fn test_ty() {
        assert_eq!(ty(&Key::ReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::ForwardStrandReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::ReverseStrandReadDepths), Some(Type::Integer));
        assert_eq!(ty(&Key::ReadDepth), Some(Type::Integer));
        assert_eq!(ty(&Key::ExpectedAlternateAlleleCounts), Some(Type::Integer));
        assert_eq!(ty(&Key::Filter), Some(Type::String));
        assert_eq!(ty(&Key::GenotypeLikelihoods), Some(Type::Float));
        assert_eq!(ty(&Key::GenotypePosteriorProbabilities), Some(Type::Float));
        assert_eq!(ty(&Key::ConditionalGenotypeQuality), Some(Type::Integer));
        assert_eq!(ty(&Key::Genotype), Some(Type::String));
        assert_eq!(ty(&Key::HaplotypeQuality), Some(Type::Integer));
        assert_eq!(ty(&Key::MappingQuality), Some(Type::Integer));
        assert_eq!(ty(&Key::RoundedGenotypeLikelihoods), Some(Type::Integer));
        assert_eq!(
            ty(&Key::RoundedGenotypePosteriorProbabilities),
            Some(Type::Integer)
        );
        assert_eq!(ty(&Key::PhasingQuality), Some(Type::Integer));
        assert_eq!(ty(&Key::PhaseSet), Some(Type::Integer));

        assert_eq!(ty(&Key::GenotypeCopyNumber), Some(Type::Integer));
        assert_eq!(ty(&Key::GenotypeCopyNumberQuality), Some(Type::Float));
        assert_eq!(ty(&Key::GenotypeCopyNumberLikelihoods), Some(Type::Float));
        assert_eq!(
            ty(&Key::GenotypeCopyNumberPosteriorProbabilities),
            Some(Type::Float)
        );
        assert_eq!(ty(&Key::NovelVariantQualityScore), Some(Type::Integer));
        assert_eq!(ty(&Key::HaplotypeId), Some(Type::Integer));
        assert_eq!(ty(&Key::AncestralHaplotypeId), Some(Type::Integer));

        assert!(ty(&Key::Other(String::from("NDLS"))).is_none());
    }

    #[test]
    fn test_description() {
        assert_eq!(
            description(&Key::ReadDepths),
            Some("Read depth for each allele")
        );
        assert_eq!(
            description(&Key::ForwardStrandReadDepths),
            Some("Read depth for each allele on the forward strand")
        );
        assert_eq!(
            description(&Key::ReverseStrandReadDepths),
            Some("Read depth for each allele on the reverse strand")
        );
        assert_eq!(description(&Key::ReadDepth), Some("Read depth"));
        assert_eq!(
            description(&Key::ExpectedAlternateAlleleCounts),
            Some("Expected alternate allele counts")
        );
        assert_eq!(
            description(&Key::Filter),
            Some(r#"Filter indicating if this genotype was "called""#)
        );
        assert_eq!(
            description(&Key::GenotypeLikelihoods),
            Some("Genotype likelihoods")
        );
        assert_eq!(
            description(&Key::GenotypePosteriorProbabilities),
            Some("Genotype posterior probabilities")
        );
        assert_eq!(
            description(&Key::ConditionalGenotypeQuality),
            Some("Conditional genotype quality")
        );
        assert_eq!(description(&Key::Genotype), Some("Genotype"));
        assert_eq!(
            description(&Key::HaplotypeQuality),
            Some("Haplotype quality")
        );
        assert_eq!(
            description(&Key::MappingQuality),
            Some("RMS mapping quality")
        );
        assert_eq!(
            description(&Key::RoundedGenotypeLikelihoods),
            Some("Phred-scaled genotype likelihoods rounded to the closest integer")
        );
        assert_eq!(
            description(&Key::RoundedGenotypePosteriorProbabilities),
            Some("Phred-scaled genotype posterior probabilities rounded to the closest integer")
        );
        assert_eq!(description(&Key::PhasingQuality), Some("Phasing quality"));
        assert_eq!(description(&Key::PhaseSet), Some("Phase set"));

        assert_eq!(
            description(&Key::GenotypeCopyNumber),
            Some("Copy number genotype for imprecise events")
        );
        assert_eq!(
            description(&Key::GenotypeCopyNumberQuality),
            Some("Copy number genotype quality for imprecise events")
        );
        assert_eq!(
            description(&Key::GenotypeCopyNumberLikelihoods),
            Some("Copy number genotype likelihood for imprecise events")
        );
        assert_eq!(
            description(&Key::GenotypeCopyNumberPosteriorProbabilities),
            Some("Copy number posterior probabilities")
        );
        assert_eq!(
            description(&Key::NovelVariantQualityScore),
            Some("Phred style probability score that the variant is novel")
        );
        assert_eq!(
            description(&Key::HaplotypeId),
            Some("Unique haplotype identifier")
        );
        assert_eq!(
            description(&Key::AncestralHaplotypeId),
            Some("Unique identifier of ancestral haplotype")
        );

        assert!(description(&Key::Other(String::from("NDLS"))).is_none());
    }
}
