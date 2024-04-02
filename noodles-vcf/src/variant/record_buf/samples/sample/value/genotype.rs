//! Variant record samples genotype value.

pub mod allele;
mod parser;

pub use self::{allele::Allele, parser::ParseError};
use crate::variant::record::samples::series::value::genotype::Phasing;

use std::{io, str::FromStr};

/// A variant record samples genotype value.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Genotype(Vec<Allele>);

impl AsRef<[Allele]> for Genotype {
    fn as_ref(&self) -> &[Allele] {
        &self.0
    }
}

impl AsMut<Vec<Allele>> for Genotype {
    fn as_mut(&mut self) -> &mut Vec<Allele> {
        &mut self.0
    }
}

impl FromStr for Genotype {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parser::parse(s)
    }
}

impl Extend<Allele> for Genotype {
    fn extend<T: IntoIterator<Item = Allele>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<Allele> for Genotype {
    fn from_iter<T: IntoIterator<Item = Allele>>(iter: T) -> Self {
        let mut genotype = Self::default();
        genotype.extend(iter);
        genotype
    }
}

impl TryFrom<&dyn crate::variant::record::samples::series::value::Genotype> for Genotype {
    type Error = io::Error;

    fn try_from(
        genotype: &dyn crate::variant::record::samples::series::value::Genotype,
    ) -> Result<Self, Self::Error> {
        genotype
            .iter()
            .map(|result| result.map(|(position, phasing)| Allele::new(position, phasing)))
            .collect::<io::Result<_>>()
            .map(Self)
    }
}

impl crate::variant::record::samples::series::value::Genotype for &Genotype {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, Phasing)>> + '_> {
        Box::new(
            self.0
                .iter()
                .map(|allele| Ok((allele.position(), allele.phasing()))),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
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
}
