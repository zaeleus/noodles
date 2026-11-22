use std::{error, fmt};

use super::MISSING;
use crate::variant::record_buf::AlternateBases;

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

    for t in s.split(DELIMITER) {
        let allele = match t {
            MISSING => None,
            _ => Some(t.into()),
        };

        alleles.push(allele);
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
            [Some(String::from("A"))].into_iter().collect()
        );

        alternate_bases.as_mut().clear();
        parse_alternate_bases(".", &mut alternate_bases)?;
        assert_eq!(alternate_bases, [None].into_iter().collect());

        alternate_bases.as_mut().clear();
        parse_alternate_bases("A,C", &mut alternate_bases)?;
        assert_eq!(
            alternate_bases,
            [Some(String::from("A")), Some(String::from("C"))]
                .into_iter()
                .collect()
        );

        alternate_bases.as_mut().clear();
        parse_alternate_bases("A,.", &mut alternate_bases)?;
        assert_eq!(
            alternate_bases,
            [Some(String::from("A")), None].into_iter().collect()
        );

        alternate_bases.as_mut().clear();
        assert_eq!(
            parse_alternate_bases("", &mut alternate_bases),
            Err(ParseError::Empty)
        );

        Ok(())
    }
}
