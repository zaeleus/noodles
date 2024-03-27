use std::{error, fmt};

use super::{allele, Allele, Genotype};
use crate::variant::record::samples::series::value::genotype::Phasing;

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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FirstPhasing {
    Explicit(Phasing),
    Implicit(Phasing),
}

impl FirstPhasing {
    fn phasing(&self) -> Phasing {
        match self {
            FirstPhasing::Explicit(phasing) => *phasing,
            FirstPhasing::Implicit(phasing) => *phasing,
        }
    }
}

pub(super) fn parse(mut s: &str) -> Result<Genotype, ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    let raw_allele = next_allele(&mut s);
    let first_allele = parse_first_allele(raw_allele).map_err(ParseError::InvalidAllele)?;
    let (first_position, mut first_phasing) = match first_allele {
        (position, Some(phasing)) => (position, FirstPhasing::Explicit(phasing)),
        (position, None) => (position, FirstPhasing::Implicit(Phasing::Phased)),
    };

    let mut alleles = vec![Allele::new(first_position, first_phasing.phasing())];

    while !s.is_empty() {
        let raw_allele = next_allele(&mut s);
        let allele: Allele = raw_allele.parse().map_err(ParseError::InvalidAllele)?;

        if first_phasing == FirstPhasing::Implicit(Phasing::Phased)
            && allele.phasing() == Phasing::Unphased
        {
            first_phasing = FirstPhasing::Implicit(Phasing::Unphased);
        }

        alleles.push(allele);
    }

    *alleles[0].phasing_mut() = first_phasing.phasing();

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

fn parse_first_allele(s: &str) -> Result<(Option<usize>, Option<Phasing>), allele::ParseError> {
    use super::allele::{parse_phasing, parse_position};

    match parse_phasing(&s[..1]) {
        Ok(phasing) => {
            let position = parse_position(&s[1..])?;
            Ok((position, Some(phasing)))
        }
        Err(_) => {
            if let Ok(position) = parse_position(s) {
                Ok((position, None))
            } else {
                Err(allele::ParseError::InvalidPhasing)
            }
        }
    }
}
