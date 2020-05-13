use crate::header::{info::Type, Number};

use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Key {
    AncestralAllele,
    AlleleCount,
    TotalReadDepths,
    ForwardStrandReadDepths,
    ReverseStrandReadDepths,
    AlleleFrequencies,
    TotalAlleleCount,
    BaseQuality,
    Cigar,
    IsInDbSnp,
    TotalDepth,
    EndPosition,
    IsInHapMap2,
    IsInHapMap3,
    MappingQuality,
    ZeroMappingQualityCount,
    SamplesWithDataCount,
    StrandBias,
    IsSomaticMutation,
    IsValidated,
    IsIn1000Genomes,
    Other(String, Number, Type),
}

impl Key {
    pub fn number(&self) -> Number {
        match self {
            Self::AncestralAllele => Number::Count(1),
            Self::AlleleCount => Number::A,
            Self::TotalReadDepths => Number::R,
            Self::ForwardStrandReadDepths => Number::R,
            Self::ReverseStrandReadDepths => Number::R,
            Self::AlleleFrequencies => Number::A,
            Self::TotalAlleleCount => Number::Count(1),
            Self::BaseQuality => Number::Count(1),
            Self::Cigar => Number::A,
            Self::IsInDbSnp => Number::Count(0),
            Self::TotalDepth => Number::Count(1),
            Self::EndPosition => Number::Count(1),
            Self::IsInHapMap2 => Number::Count(0),
            Self::IsInHapMap3 => Number::Count(0),
            Self::MappingQuality => Number::Count(1),
            Self::ZeroMappingQualityCount => Number::Count(1),
            Self::SamplesWithDataCount => Number::Count(1),
            Self::StrandBias => Number::Count(4),
            Self::IsSomaticMutation => Number::Count(0),
            Self::IsValidated => Number::Count(0),
            Self::IsIn1000Genomes => Number::Count(0),
            Self::Other(_, number, _) => *number,
        }
    }

    pub fn ty(&self) -> Type {
        match self {
            Self::AncestralAllele => Type::String,
            Self::AlleleCount => Type::Integer,
            Self::TotalReadDepths => Type::Integer,
            Self::ForwardStrandReadDepths => Type::Integer,
            Self::ReverseStrandReadDepths => Type::Integer,
            Self::AlleleFrequencies => Type::Float,
            Self::TotalAlleleCount => Type::Integer,
            Self::BaseQuality => Type::Float,
            Self::Cigar => Type::String,
            Self::IsInDbSnp => Type::Flag,
            Self::TotalDepth => Type::Integer,
            Self::EndPosition => Type::Integer,
            Self::IsInHapMap2 => Type::Flag,
            Self::IsInHapMap3 => Type::Flag,
            Self::MappingQuality => Type::Float,
            Self::ZeroMappingQualityCount => Type::Integer,
            Self::SamplesWithDataCount => Type::Integer,
            Self::StrandBias => Type::Integer,
            Self::IsSomaticMutation => Type::Flag,
            Self::IsValidated => Type::Flag,
            Self::IsIn1000Genomes => Type::Flag,
            Self::Other(_, _, ty) => *ty,
        }
    }
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::AncestralAllele => "AA",
            Self::AlleleCount => "AC",
            Self::TotalReadDepths => "AD",
            Self::ForwardStrandReadDepths => "ADF",
            Self::ReverseStrandReadDepths => "ADR",
            Self::AlleleFrequencies => "AF",
            Self::TotalAlleleCount => "AN",
            Self::BaseQuality => "BQ",
            Self::Cigar => "CIGAR",
            Self::IsInDbSnp => "DB",
            Self::TotalDepth => "DP",
            Self::EndPosition => "END",
            Self::IsInHapMap2 => "H2",
            Self::IsInHapMap3 => "H3",
            Self::MappingQuality => "MQ",
            Self::ZeroMappingQualityCount => "MQ0",
            Self::SamplesWithDataCount => "NS",
            Self::StrandBias => "SB",
            Self::IsSomaticMutation => "SOMATIC",
            Self::IsValidated => "VALIDATED",
            Self::IsIn1000Genomes => "1000G",
            Self::Other(key, ..) => key,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid info key: {}", self.0)
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError(s.into()));
        }

        match s {
            "AA" => Ok(Self::AncestralAllele),
            "AC" => Ok(Self::AlleleCount),
            "AD" => Ok(Self::TotalReadDepths),
            "ADF" => Ok(Self::ForwardStrandReadDepths),
            "ADR" => Ok(Self::ReverseStrandReadDepths),
            "AF" => Ok(Self::AlleleFrequencies),
            "AN" => Ok(Self::TotalAlleleCount),
            "BQ" => Ok(Self::BaseQuality),
            "CIGAR" => Ok(Self::Cigar),
            "DB" => Ok(Self::IsInDbSnp),
            "DP" => Ok(Self::TotalDepth),
            "END" => Ok(Self::EndPosition),
            "H2" => Ok(Self::IsInHapMap2),
            "H3" => Ok(Self::IsInHapMap3),
            "MQ" => Ok(Self::MappingQuality),
            "MQ0" => Ok(Self::ZeroMappingQualityCount),
            "NS" => Ok(Self::SamplesWithDataCount),
            "SB" => Ok(Self::StrandBias),
            "SOMATIC" => Ok(Self::IsSomaticMutation),
            "VALIDATED" => Ok(Self::IsValidated),
            "1000G" => Ok(Self::IsIn1000Genomes),
            _ => Ok(Self::Other(s.into(), Number::Count(1), Type::String)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_number() {
        assert_eq!(Key::AncestralAllele.number(), Number::Count(1));
        assert_eq!(Key::AlleleCount.number(), Number::A);
        assert_eq!(Key::TotalReadDepths.number(), Number::R);
        assert_eq!(Key::ForwardStrandReadDepths.number(), Number::R);
        assert_eq!(Key::ReverseStrandReadDepths.number(), Number::R);
        assert_eq!(Key::AlleleFrequencies.number(), Number::A);
        assert_eq!(Key::TotalAlleleCount.number(), Number::Count(1));
        assert_eq!(Key::BaseQuality.number(), Number::Count(1));
        assert_eq!(Key::Cigar.number(), Number::A);
        assert_eq!(Key::IsInDbSnp.number(), Number::Count(0));
        assert_eq!(Key::TotalDepth.number(), Number::Count(1));
        assert_eq!(Key::EndPosition.number(), Number::Count(1));
        assert_eq!(Key::IsInHapMap2.number(), Number::Count(0));
        assert_eq!(Key::IsInHapMap3.number(), Number::Count(0));
        assert_eq!(Key::MappingQuality.number(), Number::Count(1));
        assert_eq!(Key::ZeroMappingQualityCount.number(), Number::Count(1));
        assert_eq!(Key::SamplesWithDataCount.number(), Number::Count(1));
        assert_eq!(Key::StrandBias.number(), Number::Count(4));
        assert_eq!(Key::IsSomaticMutation.number(), Number::Count(0));
        assert_eq!(Key::IsValidated.number(), Number::Count(0));
        assert_eq!(Key::IsIn1000Genomes.number(), Number::Count(0));
        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).number(),
            Number::Count(1)
        );
    }

    #[test]
    fn test_ty() {
        assert_eq!(Key::AncestralAllele.ty(), Type::String);
        assert_eq!(Key::AlleleCount.ty(), Type::Integer);
        assert_eq!(Key::TotalReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ForwardStrandReadDepths.ty(), Type::Integer);
        assert_eq!(Key::ReverseStrandReadDepths.ty(), Type::Integer);
        assert_eq!(Key::AlleleFrequencies.ty(), Type::Float);
        assert_eq!(Key::TotalAlleleCount.ty(), Type::Integer);
        assert_eq!(Key::BaseQuality.ty(), Type::Float);
        assert_eq!(Key::Cigar.ty(), Type::String);
        assert_eq!(Key::IsInDbSnp.ty(), Type::Flag);
        assert_eq!(Key::TotalDepth.ty(), Type::Integer);
        assert_eq!(Key::EndPosition.ty(), Type::Integer);
        assert_eq!(Key::IsInHapMap2.ty(), Type::Flag);
        assert_eq!(Key::IsInHapMap3.ty(), Type::Flag);
        assert_eq!(Key::MappingQuality.ty(), Type::Float);
        assert_eq!(Key::ZeroMappingQualityCount.ty(), Type::Integer);
        assert_eq!(Key::SamplesWithDataCount.ty(), Type::Integer);
        assert_eq!(Key::StrandBias.ty(), Type::Integer);
        assert_eq!(Key::IsSomaticMutation.ty(), Type::Flag);
        assert_eq!(Key::IsValidated.ty(), Type::Flag);
        assert_eq!(Key::IsIn1000Genomes.ty(), Type::Flag);
        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).ty(),
            Type::String
        );
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Key::AncestralAllele.to_string(), "AA");
        assert_eq!(Key::AlleleCount.to_string(), "AC");
        assert_eq!(Key::TotalReadDepths.to_string(), "AD");
        assert_eq!(Key::ForwardStrandReadDepths.to_string(), "ADF");
        assert_eq!(Key::ReverseStrandReadDepths.to_string(), "ADR");
        assert_eq!(Key::AlleleFrequencies.to_string(), "AF");
        assert_eq!(Key::TotalAlleleCount.to_string(), "AN");
        assert_eq!(Key::BaseQuality.to_string(), "BQ");
        assert_eq!(Key::Cigar.to_string(), "CIGAR");
        assert_eq!(Key::IsInDbSnp.to_string(), "DB");
        assert_eq!(Key::TotalDepth.to_string(), "DP");
        assert_eq!(Key::EndPosition.to_string(), "END");
        assert_eq!(Key::IsInHapMap2.to_string(), "H2");
        assert_eq!(Key::IsInHapMap3.to_string(), "H3");
        assert_eq!(Key::MappingQuality.to_string(), "MQ");
        assert_eq!(Key::ZeroMappingQualityCount.to_string(), "MQ0");
        assert_eq!(Key::SamplesWithDataCount.to_string(), "NS");
        assert_eq!(Key::StrandBias.to_string(), "SB");
        assert_eq!(Key::IsSomaticMutation.to_string(), "SOMATIC");
        assert_eq!(Key::IsValidated.to_string(), "VALIDATED");
        assert_eq!(Key::IsIn1000Genomes.to_string(), "1000G");
        assert_eq!(
            Key::Other(String::from("NDLS"), Number::Count(1), Type::String).to_string(),
            "NDLS"
        );
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("AA".parse::<Key>()?, Key::AncestralAllele);
        assert_eq!("AC".parse::<Key>()?, Key::AlleleCount);
        assert_eq!("AD".parse::<Key>()?, Key::TotalReadDepths);
        assert_eq!("ADF".parse::<Key>()?, Key::ForwardStrandReadDepths);
        assert_eq!("ADR".parse::<Key>()?, Key::ReverseStrandReadDepths);
        assert_eq!("AF".parse::<Key>()?, Key::AlleleFrequencies);
        assert_eq!("AN".parse::<Key>()?, Key::TotalAlleleCount);
        assert_eq!("BQ".parse::<Key>()?, Key::BaseQuality);
        assert_eq!("CIGAR".parse::<Key>()?, Key::Cigar);
        assert_eq!("DB".parse::<Key>()?, Key::IsInDbSnp);
        assert_eq!("DP".parse::<Key>()?, Key::TotalDepth);
        assert_eq!("END".parse::<Key>()?, Key::EndPosition);
        assert_eq!("H2".parse::<Key>()?, Key::IsInHapMap2);
        assert_eq!("H3".parse::<Key>()?, Key::IsInHapMap3);
        assert_eq!("MQ".parse::<Key>()?, Key::MappingQuality);
        assert_eq!("MQ0".parse::<Key>()?, Key::ZeroMappingQualityCount);
        assert_eq!("NS".parse::<Key>()?, Key::SamplesWithDataCount);
        assert_eq!("SB".parse::<Key>()?, Key::StrandBias);
        assert_eq!("SOMATIC".parse::<Key>()?, Key::IsSomaticMutation);
        assert_eq!("VALIDATED".parse::<Key>()?, Key::IsValidated);
        assert_eq!("1000G".parse::<Key>()?, Key::IsIn1000Genomes);

        assert!("".parse::<Key>().is_err());

        Ok(())
    }
}
