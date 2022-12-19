//! VCF record alternate bases allele and symbol.

pub mod symbol;

pub use self::symbol::Symbol;

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

use crate::record::reference_bases::{base, Base};

/// A VCF alternate bases allele.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Allele {
    /// A list of bases (e.g., `A`, `AC`, etc.).
    Bases(Vec<Base>),
    /// A symbolic allele (e.g., `<DEL>`, `<CN:0>`, etc.).
    Symbol(Symbol),
    /// A breakend (e.g., `]sq0:5]A`, `G.`, etc.).
    Breakend(String),
    /// An overlapping deletion, i.e., a missing allele (`*`).
    OverlappingDeletion,
}

impl fmt::Display for Allele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Bases(bases) => {
                for base in bases {
                    f.write_char(char::from(*base))?;
                }

                Ok(())
            }
            Self::Symbol(symbol) => write!(f, "<{}>", symbol),
            Self::Breakend(breakend) => f.write_str(breakend),
            Self::OverlappingDeletion => f.write_str("*"),
        }
    }
}

/// An error returned when a raw alternate bases allele fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The symbol is invalid.
    InvalidSymbol(symbol::ParseError),
    /// A base is invalid.
    InvalidBase(base::TryFromCharError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::InvalidSymbol(e) => Some(e),
            Self::InvalidBase(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidSymbol(_) => f.write_str("invalid symbol"),
            Self::InvalidBase(_) => f.write_str("invalid base"),
        }
    }
}

impl FromStr for Allele {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "*" => Ok(Self::OverlappingDeletion),
            _ => {
                if s.starts_with('<') {
                    s.trim_matches(|c| c == '<' || c == '>')
                        .parse()
                        .map(Self::Symbol)
                        .map_err(ParseError::InvalidSymbol)
                } else if is_breakend(s) {
                    Ok(Self::Breakend(s.into()))
                } else {
                    s.chars()
                        .map(|c| c.to_ascii_uppercase())
                        .map(Base::try_from)
                        .collect::<Result<_, _>>()
                        .map(Self::Bases)
                        .map_err(ParseError::InvalidBase)
                }
            }
        }
    }
}

fn is_breakend(s: &str) -> bool {
    s.contains(|c| c == '[' || c == ']') || s.starts_with('.') || s.ends_with('.')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let allele = Allele::Bases(vec![Base::G]);
        assert_eq!(allele.to_string(), "G");

        let allele = Allele::Bases(vec![Base::G, Base::T]);
        assert_eq!(allele.to_string(), "GT");

        let allele = Allele::Symbol(Symbol::StructuralVariant(symbol::StructuralVariant::from(
            symbol::structural_variant::Type::Duplication,
        )));
        assert_eq!(allele.to_string(), "<DUP>");

        let allele = Allele::Symbol(Symbol::NonstructuralVariant(String::from("CN0")));
        assert_eq!(allele.to_string(), "<CN0>");

        let allele = Allele::Symbol(Symbol::NonstructuralVariant(String::from("CN:0")));
        assert_eq!(allele.to_string(), "<CN:0>");

        let allele = Allele::Breakend(String::from("]sq0:5]A"));
        assert_eq!(allele.to_string(), "]sq0:5]A");

        let allele = Allele::Breakend(String::from("C[sq1:13["));
        assert_eq!(allele.to_string(), "C[sq1:13[");

        let allele = Allele::Breakend(String::from("G."));
        assert_eq!(allele.to_string(), "G.");

        let allele = Allele::Breakend(String::from("CT."));
        assert_eq!(allele.to_string(), "CT.");

        let allele = Allele::Breakend(String::from(".A"));
        assert_eq!(allele.to_string(), ".A");

        let allele = Allele::Breakend(String::from(".GC"));
        assert_eq!(allele.to_string(), ".GC");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("G".parse::<Allele>(), Ok(Allele::Bases(vec![Base::G])));

        assert_eq!(
            "GT".parse::<Allele>(),
            Ok(Allele::Bases(vec![Base::G, Base::T]))
        );

        assert_eq!(
            "<DUP>".parse::<Allele>(),
            Ok(Allele::Symbol(Symbol::StructuralVariant(
                symbol::StructuralVariant::from(symbol::structural_variant::Type::Duplication,)
            )))
        );

        assert_eq!(
            "<CN0>".parse::<Allele>(),
            Ok(Allele::Symbol(Symbol::NonstructuralVariant(String::from(
                "CN0"
            ))))
        );

        assert_eq!(
            "<CN:0>".parse::<Allele>(),
            Ok(Allele::Symbol(Symbol::NonstructuralVariant(String::from(
                "CN:0"
            ))))
        );

        assert_eq!(
            "]sq0:5]A".parse::<Allele>(),
            Ok(Allele::Breakend(String::from("]sq0:5]A")))
        );

        assert_eq!(
            "C[sq1:13[".parse::<Allele>(),
            Ok(Allele::Breakend(String::from("C[sq1:13[")))
        );

        assert_eq!(
            "G.".parse::<Allele>(),
            Ok(Allele::Breakend(String::from("G.")))
        );

        assert_eq!(
            "CT.".parse::<Allele>(),
            Ok(Allele::Breakend(String::from("CT.")))
        );

        assert_eq!(
            ".A".parse::<Allele>(),
            Ok(Allele::Breakend(String::from(".A")))
        );

        assert_eq!(
            ".GC".parse::<Allele>(),
            Ok(Allele::Breakend(String::from(".GC")))
        );

        assert_eq!("".parse::<Allele>(), Err(ParseError::Empty));
        assert!(matches!(
            "<>".parse::<Allele>(),
            Err(ParseError::InvalidSymbol(_))
        ));
        assert!(matches!(
            "Z".parse::<Allele>(),
            Err(ParseError::InvalidBase(_))
        ));
    }
}
