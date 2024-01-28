use std::{error, fmt};

use crate::record::AlternateBases;

/// An error when raw VCF record alternate bases fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
        }
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

    let alleles = alternate_bases.as_mut();

    for allele in s.split(DELIMITER) {
        alleles.push(allele.into());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_alternate_bases() -> Result<(), ParseError> {
        let mut alternate_bases = AlternateBases::default();

        alternate_bases.as_mut().clear();
        parse_alternate_bases("A", &mut alternate_bases)?;
        assert_eq!(
            alternate_bases,
            AlternateBases::from(vec![String::from("A")])
        );

        alternate_bases.as_mut().clear();
        parse_alternate_bases("A,C", &mut alternate_bases)?;
        assert_eq!(
            alternate_bases,
            AlternateBases::from(vec![String::from("A"), String::from("C")])
        );

        alternate_bases.as_mut().clear();
        assert_eq!(
            parse_alternate_bases("", &mut alternate_bases),
            Err(ParseError::Empty)
        );

        Ok(())
    }
}
