use std::{error, fmt};

use super::{allele, Genotype};

/// An error returned when a raw VCF record genotype value fails to parse.
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
            Self::Empty => f.write_str("empty input"),
            Self::InvalidAllele(_) => f.write_str("invalid allele"),
        }
    }
}

pub(super) fn parse(mut s: &str) -> Result<Genotype, ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    let mut alleles = Vec::new();

    while !s.is_empty() {
        let raw_allele = next_allele(&mut s);
        let allele = raw_allele.parse().map_err(ParseError::InvalidAllele)?;
        alleles.push(allele);
    }

    Ok(Genotype(alleles))
}

fn next_allele<'a>(s: &mut &'a str) -> &'a str {
    let (t, rest) = match s.chars().skip(1).position(is_phasing_indicator) {
        Some(i) => s.split_at(i + 1),
        None => s.split_at(s.len()),
    };

    *s = rest;

    t
}

fn is_phasing_indicator(c: char) -> bool {
    matches!(c, '/' | '|')
}
