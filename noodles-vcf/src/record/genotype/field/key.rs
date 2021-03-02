//! VCF record genotype field key.

use std::{error, fmt, str::FromStr};

use crate::header::{format::Type, Number};

/// A VCF record genotype field key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Key {
    // § 1.6.2 Genotype fields (2020-04-02)
    /// (`AD`).
    ReadDepths,
    /// (`ADF`).
    ForwardStrandReadDepths,
    /// (`ADR`).
    ReverseStrandReadDepths,
    /// (`DP`).
    ReadDepth,
    /// (`EC`).
    ExpectedAlternateAlleleCounts,
    /// (`FT`).
    Filter,
    /// (`GL`).
    GenotypeLikelihoods,
    /// (`GP`).
    GenotypePosteriorProbabilities,
    /// (`GQ`).
    ConditionalGenotypeQuality,
    /// (`GT`).
    Genotype,
    /// (`HQ`).
    HaplotypeQuality,
    /// (`MQ`).
    MappingQuality,
    /// (`PL`).
    RoundedGenotypeLikelihoods,
    /// (`PP`).
    RoundedGenotypePosteriorProbabilities,
    /// (`PQ`).
    PhasingQuality,
    /// (`PS`).
    PhaseSet,

    // § 4 FORMAT keys used for structural variants (2020-04-02)
    /// (`CN`).
    GenotypeCopyNumber,
    /// (`CNQ`).
    GenotypeCopyNumberQuality,
    /// (`CNL`).
    GenotypeCopyNumberLikelihoods,
    /// (`CNP`).
    GenotypeCopyNumberPosteriorProbabilities,
    /// (`NQ`).
    NovelVariantQualityScore,
    /// (`HAP`).
    HaplotypeId,
    /// (`AHAP`).
    AncestralHaplotypeId,

    /// Any other non-reserved key.
    Other(String, Number, Type),
}

impl Key {
    /// Returns the cardinality of the genotype field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{header::Number, record::genotype::field::Key};
    /// assert_eq!(Key::Genotype.number(), Number::Count(1));
    /// ```
    pub fn number(&self) -> Number {
        match self {
            Self::ReadDepths => Number::R,
            Self::ForwardStrandReadDepths => Number::R,
            Self::ReverseStrandReadDepths => Number::R,
            Self::ReadDepth => Number::Count(1),
            Self::ExpectedAlternateAlleleCounts => Number::A,
            Self::Filter => Number::Count(1),
            Self::GenotypeLikelihoods => Number::G,
            Self::GenotypePosteriorProbabilities => Number::G,
            Self::ConditionalGenotypeQuality => Number::Count(1),
            Self::Genotype => Number::Count(1),
            Self::HaplotypeQuality => Number::Count(2),
            Self::MappingQuality => Number::Count(1),
            Self::RoundedGenotypeLikelihoods => Number::G,
            Self::RoundedGenotypePosteriorProbabilities => Number::G,
            Self::PhasingQuality => Number::Count(1),
            Self::PhaseSet => Number::Count(1),

            Self::GenotypeCopyNumber => Number::Count(1),
            Self::GenotypeCopyNumberQuality => Number::Count(1),
            Self::GenotypeCopyNumberLikelihoods => Number::G,
            Self::GenotypeCopyNumberPosteriorProbabilities => Number::G,
            Self::NovelVariantQualityScore => Number::Count(1),
            Self::HaplotypeId => Number::Count(1),
            Self::AncestralHaplotypeId => Number::Count(1),

            Self::Other(_, number, _) => *number,
        }
    }

    /// Returns the type of the genotype field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{header::format::Type, record::genotype::field::Key};
    /// assert_eq!(Key::Genotype.ty(), Type::String);
    /// ```
    pub fn ty(&self) -> Type {
        match self {
            Self::ReadDepths => Type::Integer,
            Self::ForwardStrandReadDepths => Type::Integer,
            Self::ReverseStrandReadDepths => Type::Integer,
            Self::ReadDepth => Type::Integer,
            Self::ExpectedAlternateAlleleCounts => Type::Integer,
            Self::Filter => Type::String,
            Self::GenotypeLikelihoods => Type::Float,
            Self::GenotypePosteriorProbabilities => Type::Float,
            Self::ConditionalGenotypeQuality => Type::Integer,
            Self::Genotype => Type::String,
            Self::HaplotypeQuality => Type::Integer,
            Self::MappingQuality => Type::Integer,
            Self::RoundedGenotypeLikelihoods => Type::Integer,
            Self::RoundedGenotypePosteriorProbabilities => Type::Integer,
            Self::PhasingQuality => Type::Integer,
            Self::PhaseSet => Type::Integer,

            Self::GenotypeCopyNumber => Type::Integer,
            Self::GenotypeCopyNumberQuality => Type::Float,
            Self::GenotypeCopyNumberLikelihoods => Type::Float,
            Self::GenotypeCopyNumberPosteriorProbabilities => Type::Float,
            Self::NovelVariantQualityScore => Type::Integer,
            Self::HaplotypeId => Type::Integer,
            Self::AncestralHaplotypeId => Type::Integer,

            Self::Other(_, _, ty) => *ty,
        }
    }
}

impl AsRef<str> for Key {
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

            Self::GenotypeCopyNumber => "CN",
            Self::GenotypeCopyNumberQuality => "CNQ",
            Self::GenotypeCopyNumberLikelihoods => "CNL",
            Self::GenotypeCopyNumberPosteriorProbabilities => "CNP",
            Self::NovelVariantQualityScore => "NQ",
            Self::HaplotypeId => "HAP",
            Self::AncestralHaplotypeId => "AHAP",

            Self::Other(key, ..) => key,
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
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
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

            "CN" => Ok(Self::GenotypeCopyNumber),
            "CNQ" => Ok(Self::GenotypeCopyNumberQuality),
            "CNL" => Ok(Self::GenotypeCopyNumberLikelihoods),
            "CNP" => Ok(Self::GenotypeCopyNumberPosteriorProbabilities),
            "NQ" => Ok(Self::NovelVariantQualityScore),
            "HAP" => Ok(Self::HaplotypeId),
            "AHAP" => Ok(Self::AncestralHaplotypeId),

            _ => Ok(Self::Other(s.into(), Number::Count(1), Type::String)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_number() {
        assert_eq!(Key::ReadDepths.number(), Number::R);
        assert_eq!(Key::ForwardStrandReadDepths.number(), Number::R);
        assert_eq!(Key::ReverseStrandReadDepths.number(), Number::R);
        assert_eq!(Key::ReadDepth.number(), Number::Count(1));
        assert_eq!(Key::ExpectedAlternateAlleleCounts.number(), Number::A);
        assert_eq!(Key::Filter.number(), Number::Count(1));
        assert_eq!(Key::GenotypeLikelihoods.number(), Number::G);
        assert_eq!(Key::GenotypePosteriorProbabilities.number(), Number::G);
        assert_eq!(Key::ConditionalGenotypeQuality.number(), Number::Count(1));
        assert_eq!(Key::Genotype.number(), Number::Count(1));
        assert_eq!(Key::HaplotypeQuality.number(), Number::Count(2));
        assert_eq!(Key::MappingQuality.number(), Number::Count(1));
        assert_eq!(Key::RoundedGenotypeLikelihoods.number(), Number::G);
        assert_eq!(
            Key::RoundedGenotypePosteriorProbabilities.number(),
            Number::G
        );
        assert_eq!(Key::PhasingQuality.number(), Number::Count(1));
        assert_eq!(Key::PhaseSet.number(), Number::Count(1));

        assert_eq!(Key::GenotypeCopyNumber.number(), Number::Count(1));
        assert_eq!(Key::GenotypeCopyNumberQuality.number(), Number::Count(1));
        assert_eq!(Key::GenotypeCopyNumberLikelihoods.number(), Number::G);
        assert_eq!(
            Key::GenotypeCopyNumberPosteriorProbabilities.number(),
            Number::G
        );
        assert_eq!(Key::NovelVariantQualityScore.number(), Number::Count(1));
        assert_eq!(Key::HaplotypeId.number(), Number::Count(1));
        assert_eq!(Key::AncestralHaplotypeId.number(), Number::Count(1));

        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).number(),
            Number::Count(1)
        );
    }

    #[test]
    fn test_ty() {
        assert_eq!(Key::ReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ForwardStrandReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ReverseStrandReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ReadDepth.ty(), Type::Integer);
        assert_eq!(Key::ExpectedAlternateAlleleCounts.ty(), Type::Integer);
        assert_eq!(Key::Filter.ty(), Type::String);
        assert_eq!(Key::GenotypeLikelihoods.ty(), Type::Float);
        assert_eq!(Key::GenotypePosteriorProbabilities.ty(), Type::Float);
        assert_eq!(Key::ConditionalGenotypeQuality.ty(), Type::Integer);
        assert_eq!(Key::Genotype.ty(), Type::String);
        assert_eq!(Key::HaplotypeQuality.ty(), Type::Integer);
        assert_eq!(Key::MappingQuality.ty(), Type::Integer);
        assert_eq!(Key::RoundedGenotypeLikelihoods.ty(), Type::Integer);
        assert_eq!(
            Key::RoundedGenotypePosteriorProbabilities.ty(),
            Type::Integer
        );
        assert_eq!(Key::PhasingQuality.ty(), Type::Integer);
        assert_eq!(Key::PhaseSet.ty(), Type::Integer);

        assert_eq!(Key::GenotypeCopyNumber.ty(), Type::Integer);
        assert_eq!(Key::GenotypeCopyNumberQuality.ty(), Type::Float);
        assert_eq!(Key::GenotypeCopyNumberLikelihoods.ty(), Type::Float);
        assert_eq!(
            Key::GenotypeCopyNumberPosteriorProbabilities.ty(),
            Type::Float
        );
        assert_eq!(Key::NovelVariantQualityScore.ty(), Type::Integer);
        assert_eq!(Key::HaplotypeId.ty(), Type::Integer);
        assert_eq!(Key::AncestralHaplotypeId.ty(), Type::Integer);

        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).ty(),
            Type::String
        );
    }

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

        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).to_string(),
            "NDLS"
        );
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

        assert_eq!(
            "NDLS".parse(),
            Ok(Key::Other(
                String::from("NDLS"),
                Number::Count(1),
                Type::String
            ))
        );

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
    }
}
