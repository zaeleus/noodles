//! VCF record alternate bases allele structural variant symbol type.

use std::{error, fmt, str::FromStr};

/// A VCF alternate bases allele structural variant symbol type.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Type {
    /// A deletion (`DEL`).
    Deletion,
    /// A insertion (`INS`).
    Insertion,
    /// A duplication (`DUP`).
    Duplication,
    /// An inversion (`INV`).
    Inversion,
    /// A copy number variation (`CNV`).
    CopyNumberVariation,
    /// A breakend (`BND`).
    ///
    /// Removed in VCF 4.4.
    Breakend,
}

impl AsRef<str> for Type {
    fn as_ref(&self) -> &str {
        match self {
            Self::Deletion => "DEL",
            Self::Insertion => "INS",
            Self::Duplication => "DUP",
            Self::Inversion => "INV",
            Self::CopyNumberVariation => "CNV",
            Self::Breakend => "BND",
        }
    }
}

impl fmt::Display for Type {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF record a structural variant type of an alternate base allele
/// symbol fails to parse.
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

impl FromStr for Type {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "DEL" => Ok(Self::Deletion),
            "INS" => Ok(Self::Insertion),
            "DUP" => Ok(Self::Duplication),
            "INV" => Ok(Self::Inversion),
            "CNV" => Ok(Self::CopyNumberVariation),
            "BND" => Ok(Self::Breakend),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Type::Deletion.to_string(), "DEL");
        assert_eq!(Type::Insertion.to_string(), "INS");
        assert_eq!(Type::Duplication.to_string(), "DUP");
        assert_eq!(Type::Inversion.to_string(), "INV");
        assert_eq!(Type::CopyNumberVariation.to_string(), "CNV");
        assert_eq!(Type::Breakend.to_string(), "BND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("DEL".parse(), Ok(Type::Deletion));
        assert_eq!("INS".parse(), Ok(Type::Insertion));
        assert_eq!("DUP".parse(), Ok(Type::Duplication));
        assert_eq!("INV".parse(), Ok(Type::Inversion));
        assert_eq!("CNV".parse(), Ok(Type::CopyNumberVariation));
        assert_eq!("BND".parse(), Ok(Type::Breakend));

        assert_eq!("".parse::<Type>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Type>(), Err(ParseError::Invalid));
    }
}
