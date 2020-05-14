use std::{convert::TryFrom, error, fmt, str::FromStr};

use crate::record::reference_bases::Base;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Allele {
    Bases(Vec<Base>),
    Symbol(String),
    Breakend(String),
    OverlappingDeletion,
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid alternate bases allele: {}", self.0)
    }
}

impl FromStr for Allele {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            "*" => Ok(Self::OverlappingDeletion),
            _ => {
                if s.starts_with('<') {
                    Ok(Self::Symbol(s.into()))
                } else if s.contains(|c| c == '[' || c == ']') || (s.len() == 2 && s.contains('.'))
                {
                    Ok(Self::Breakend(s.into()))
                } else {
                    s.chars()
                        .map(|c| c.to_ascii_uppercase())
                        .map(Base::try_from)
                        .collect::<Result<_, _>>()
                        .map(Self::Bases)
                        .map_err(|_| ParseError(s.into()))
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("G".parse::<Allele>()?, Allele::Bases(vec![Base::G]));

        assert_eq!(
            "GT".parse::<Allele>()?,
            Allele::Bases(vec![Base::G, Base::T])
        );

        assert_eq!(
            "<DUP>".parse::<Allele>()?,
            Allele::Symbol(String::from("<DUP>"))
        );

        assert_eq!(
            "]sq0:5]A".parse::<Allele>()?,
            Allele::Breakend(String::from("]sq0:5]A"))
        );

        assert_eq!(
            "C[sq1:13[".parse::<Allele>()?,
            Allele::Breakend(String::from("C[sq1:13["))
        );

        assert_eq!(
            "G.".parse::<Allele>()?,
            Allele::Breakend(String::from("G."))
        );

        assert_eq!(
            ".A".parse::<Allele>()?,
            Allele::Breakend(String::from(".A"))
        );

        assert!("".parse::<Allele>().is_err());

        Ok(())
    }
}
