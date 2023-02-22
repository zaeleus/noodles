//! VCF record genotype value.

pub mod allele;
mod parser;

pub use self::{allele::Allele, parser::ParseError};

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

impl FromStr for Genotype {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parser::parse(s)
    }
}

impl TryFrom<Vec<Allele>> for Genotype {
    type Error = TryFromAllelesError;

    fn try_from(alleles: Vec<Allele>) -> Result<Self, Self::Error> {
        if alleles.is_empty() {
            Err(TryFromAllelesError::Empty)
        } else {
            Ok(Self(alleles))
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
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
            ]))
        );

        assert_eq!(
            "0|1".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), Phasing::Phased),
                Allele::new(Some(1), Phasing::Phased),
            ]))
        );

        assert_eq!(
            "./.".parse(),
            Ok(Genotype(vec![
                Allele::new(None, Phasing::Unphased),
                Allele::new(None, Phasing::Unphased),
            ]))
        );

        assert_eq!(
            "0".parse(),
            Ok(Genotype(vec![Allele::new(Some(0), Phasing::Phased)]))
        );

        assert_eq!(
            "0/1/2".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Unphased),
            ]))
        );

        assert_eq!(
            "0/1|2".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), Phasing::Unphased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Phased),
            ]))
        );

        assert_eq!(
            "|0/1/2".parse(),
            Ok(Genotype(vec![
                Allele::new(Some(0), Phasing::Phased),
                Allele::new(Some(1), Phasing::Unphased),
                Allele::new(Some(2), Phasing::Unphased),
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
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
        ]);
        assert_eq!(Genotype::try_from(expected.0.clone()), Ok(expected));

        assert_eq!(
            Genotype::try_from(Vec::new()),
            Err(TryFromAllelesError::Empty)
        );
    }
}
