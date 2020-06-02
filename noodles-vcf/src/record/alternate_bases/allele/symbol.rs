pub mod structural_variant;

pub use self::structural_variant::StructuralVariant;

use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Symbol {
    StructuralVariant(StructuralVariant),
    NonstructuralVariant(String),
    Unspecified,
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid alternate bases allele symbol: {}", self.0)
    }
}

impl FromStr for Symbol {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            "*" | "NON_REF" => Ok(Self::Unspecified),
            _ => s
                .parse::<StructuralVariant>()
                .map(Self::StructuralVariant)
                .or_else(|_| Ok(Self::NonstructuralVariant(s.into()))),
        }
    }
}

impl fmt::Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::StructuralVariant(sv) => write!(f, "{}", sv),
            Self::NonstructuralVariant(nsv) => f.write_str(nsv),
            Self::Unspecified => f.write_str("*"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!(
            "DEL".parse::<Symbol>()?,
            Symbol::StructuralVariant(StructuralVariant::from(structural_variant::Type::Deletion))
        );

        assert_eq!(
            "CN:0".parse::<Symbol>()?,
            Symbol::NonstructuralVariant(String::from("CN:0"))
        );

        assert_eq!("NON_REF".parse::<Symbol>()?, Symbol::Unspecified);
        assert_eq!("*".parse::<Symbol>()?, Symbol::Unspecified);

        assert!("".parse::<Symbol>().is_err());

        Ok(())
    }

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
}
