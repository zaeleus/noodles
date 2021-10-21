//! VCF record genotype value.

pub mod allele;

pub use self::allele::Allele;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
    str::FromStr,
};

/// A VCF record genotype value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Genotype(Vec<Allele>);

impl Deref for Genotype {
    type Target = [Allele];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Genotype {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/// An error returned when a raw VCF record genotype value fails to parse.
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

impl FromStr for Genotype {
    type Err = ParseError;

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut alleles = Vec::new();

        while let Some(i) = s.chars().skip(1).position(|c| matches!(c, '/' | '|')) {
            let allele = s[..i + 1].parse().map_err(ParseError::InvalidAllele)?;
            alleles.push(allele);
            s = &s[i + 1..];
        }

        let allele = s.parse().map_err(ParseError::InvalidAllele)?;
        alleles.push(allele);

        Ok(Self(alleles))
    }
}

impl TryFrom<Vec<Allele>> for Genotype {
    type Error = TryFromAllelesError;

    fn try_from(alleles: Vec<Allele>) -> Result<Self, Self::Error> {
        match alleles.first().map(|a| a.phasing()) {
            None => Err(TryFromAllelesError::Empty),
            Some(Some(_)) => Err(TryFromAllelesError::InvalidFirstAllelePhasing),
            Some(None) => Ok(Self(alleles)),
        }
    }
}

/// An error returned when a VCF record genotype alleles fail to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromAllelesError {
    /// The list of alleles is empty.
    Empty,
    /// The phasing of the first allele is invalid.
    ///
    /// The first allele cannot have phasing information.
    InvalidFirstAllelePhasing,
}

impl error::Error for TryFromAllelesError {}

impl fmt::Display for TryFromAllelesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidFirstAllelePhasing => f.write_str("invalid first allele phasing"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        use allele::Phasing;

        assert_eq!(
            "0/1".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), None),
                Allele::new(Some(1), Some(Phasing::Unphased)),
            ]))
        );

        assert_eq!(
            "0|1".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), None),
                Allele::new(Some(1), Some(Phasing::Phased)),
            ]))
        );

        assert_eq!(
            "./.".parse(),
            Ok(Genotype(vec![
                Allele::new(None, None),
                Allele::new(None, Some(Phasing::Unphased)),
            ]))
        );

        assert_eq!("0".parse(), Ok(Genotype(vec![Allele::new(Some(0), None)])));

        assert_eq!(
            "0/1/2".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), None),
                Allele::new(Some(1), Some(Phasing::Unphased)),
                Allele::new(Some(2), Some(Phasing::Unphased)),
            ]))
        );

        assert_eq!(
            "0/1|2".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), None),
                Allele::new(Some(1), Some(Phasing::Unphased)),
                Allele::new(Some(2), Some(Phasing::Phased)),
            ]))
        );

        assert!(matches!(
            "0:1".parse::<Genotype>(),
            Err(ParseError::InvalidAllele(_))
        ));
    }

    #[test]
    fn test_try_from_alleles_for_genotype() {
        use allele::Phasing;

        let expected = Genotype(vec![
            Allele::new(Some(0), None),
            Allele::new(Some(1), Some(Phasing::Unphased)),
        ]);
        assert_eq!(Genotype::try_from(expected.0.clone()), Ok(expected));

        assert_eq!(
            Genotype::try_from(vec![Allele::new(Some(0), Some(Phasing::Unphased))]),
            Err(TryFromAllelesError::InvalidFirstAllelePhasing),
        );

        assert_eq!(
            Genotype::try_from(Vec::new()),
            Err(TryFromAllelesError::Empty)
        );
    }
}
