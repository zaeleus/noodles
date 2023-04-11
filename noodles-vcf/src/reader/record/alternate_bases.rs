use std::{error, fmt};

use noodles_core as core;

use crate::record::{alternate_bases::allele, AlternateBases};

/// An error when raw VCF record alternate bases fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// An allele is invalid.
    InvalidAllele(allele::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::InvalidAllele(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidAllele(_) => write!(f, "invalid allele"),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_alternate_bases(
    s: &str,
    alternate_bases: &mut AlternateBases,
) -> Result<(), ParseError> {
    const DELIMITER: char = ',';

    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    for raw_allele in s.split(DELIMITER) {
        let allele = raw_allele.parse().map_err(ParseError::InvalidAllele)?;
        alternate_bases.push(allele);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::{alternate_bases::Allele, reference_bases::Base};

    #[test]
    fn test_parse_alternate_bases() -> Result<(), ParseError> {
        let mut alternate_bases = AlternateBases::default();

        alternate_bases.clear();
        parse_alternate_bases("A", &mut alternate_bases)?;
        assert_eq!(
            alternate_bases,
            AlternateBases::from(vec![Allele::Bases(vec![Base::A])])
        );

        alternate_bases.clear();
        parse_alternate_bases("A,C", &mut alternate_bases)?;
        assert_eq!(
            alternate_bases,
            AlternateBases::from(vec![
                Allele::Bases(vec![Base::A]),
                Allele::Bases(vec![Base::C]),
            ])
        );

        alternate_bases.clear();
        assert_eq!(
            parse_alternate_bases("", &mut alternate_bases),
            Err(ParseError::Empty)
        );

        alternate_bases.clear();
        assert!(matches!(
            parse_alternate_bases(".", &mut alternate_bases),
            Err(ParseError::InvalidAllele(_))
        ));

        Ok(())
    }
}
