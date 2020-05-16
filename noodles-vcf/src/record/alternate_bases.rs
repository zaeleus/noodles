mod allele;

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
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid alternate bases: {}", self.0)
    }
}

impl FromStr for AlternateBases {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            MISSING_FIELD => Ok(AlternateBases::default()),
            _ => s
                .split(DELIMITER)
                .map(|s| s.parse().map_err(|_| ParseError(s.into())))
                .collect::<Result<_, _>>()
                .map(AlternateBases),
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
    fn test_from_str() -> Result<(), ParseError> {
        assert!(".".parse::<AlternateBases>()?.is_empty());

        let alternate_baes = "G".parse::<AlternateBases>()?;
        assert_eq!(alternate_baes.len(), 1);

        let alternate_baes = "G,T".parse::<AlternateBases>()?;
        assert_eq!(alternate_baes.len(), 2);

        assert!("".parse::<AlternateBases>().is_err());

        Ok(())
    }
}
