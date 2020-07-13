//! VCF record alternate bases allele symbol and structural variant.

pub mod structural_variant;

pub use self::structural_variant::StructuralVariant;

use std::{error, fmt, str::FromStr};

/// A VCF alternate bases allele symbol.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Symbol {
    /// A structural variant.
    StructuralVariant(StructuralVariant),
    /// A nonstructural variant.
    NonstructuralVariant(String),
    /// An unspecific symbol.
    Unspecified,
}

impl fmt::Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StructuralVariant(sv) => write!(f, "{}", sv),
            Self::NonstructuralVariant(nsv) => f.write_str(nsv),
            Self::Unspecified => f.write_str("*"),
        }
    }
}

/// An error returned when a raw VCF record alternate base allele symbol fails to parse.
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

impl FromStr for Symbol {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "*" | "NON_REF" => Ok(Self::Unspecified),
            _ => s
                .parse::<StructuralVariant>()
                .map(Self::StructuralVariant)
                .or_else(|_| Ok(Self::NonstructuralVariant(s.into()))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let symbol =
            Symbol::StructuralVariant(StructuralVariant::from(structural_variant::Type::Deletion));
        assert_eq!(symbol.to_string(), "DEL");

        let symbol = Symbol::NonstructuralVariant(String::from("CN:0"));
        assert_eq!(symbol.to_string(), "CN:0");

        let symbol = Symbol::Unspecified;
        assert_eq!(symbol.to_string(), "*");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "DEL".parse(),
            Ok(Symbol::StructuralVariant(StructuralVariant::from(
                structural_variant::Type::Deletion
            )))
        );

        assert_eq!(
            "CN:0".parse(),
            Ok(Symbol::NonstructuralVariant(String::from("CN:0")))
        );

        assert_eq!("NON_REF".parse(), Ok(Symbol::Unspecified));
        assert_eq!("*".parse(), Ok(Symbol::Unspecified));

        assert_eq!("".parse::<Symbol>(), Err(ParseError::Empty));
    }
}
