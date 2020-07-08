pub mod allele;

pub use self::allele::Allele;

use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

const DELIMITER: char = ',';

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlternateBases(Vec<Allele>);

impl Deref for AlternateBases {
    type Target = [Allele];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for AlternateBases {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, allele) in self.iter().enumerate() {
                if i > 0 {
                    f.write_str(",")?;
                }

                write!(f, "{}", allele)?;
            }

            Ok(())
        }
    }
}

impl From<Vec<Allele>> for AlternateBases {
    fn from(alleles: Vec<Allele>) -> Self {
        Self(alleles)
    }
}

/// An error returned when raw VCF alternate bases fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// An allele is invalid.
    InvalidAllele(allele::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidAllele(e) => write!(f, "invalid allele: {}", e),
        }
    }
}

impl FromStr for AlternateBases {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(AlternateBases::default()),
            _ => s
                .split(DELIMITER)
                .map(|s| s.parse().map_err(ParseError::InvalidAllele))
                .collect::<Result<_, _>>()
                .map(Self),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::record::reference_bases::Base;

    use super::*;

    #[test]
    fn test_fmt() {
        let alternate_bases = AlternateBases(vec![Allele::Bases(vec![Base::G])]);
        assert_eq!(alternate_bases.to_string(), "G");

        let alternate_bases = AlternateBases(vec![
            Allele::Bases(vec![Base::G]),
            Allele::Bases(vec![Base::T]),
        ]);
        assert_eq!(alternate_bases.to_string(), "G,T");

        let alternate_bases = AlternateBases(vec![]);
        assert_eq!(alternate_bases.to_string(), ".");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(".".parse(), Ok(AlternateBases::default()));

        assert_eq!(
            "G".parse(),
            Ok(AlternateBases::from(vec![Allele::Bases(vec![Base::G])]))
        );

        assert_eq!(
            "G,T".parse(),
            Ok(AlternateBases::from(vec![
                Allele::Bases(vec![Base::G]),
                Allele::Bases(vec![Base::T]),
            ]))
        );

        assert_eq!("".parse::<AlternateBases>(), Err(ParseError::Empty));
    }
}
