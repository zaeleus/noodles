//! Variant record samples genotype value allele.

use std::{error, fmt, num, str::FromStr};

use crate::variant::record::samples::series::value::genotype::Phasing;

/// A VCF record genotype value allele.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Allele {
    position: Option<usize>,
    phasing: Phasing,
}

impl Allele {
    /// Creates a variant record samples genotype value allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::series::value::genotype::Phasing,
    ///     record_buf::samples::sample::value::genotype::Allele,
    /// };
    ///
    /// let allele = Allele::new(Some(0), Phasing::Phased);
    /// ```
    pub fn new(position: Option<usize>, phasing: Phasing) -> Self {
        Self { position, phasing }
    }

    /// Returns the position of the allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::series::value::genotype::Phasing,
    ///     record_buf::samples::sample::value::genotype::Allele,
    /// };
    ///
    /// let allele = Allele::new(Some(0), Phasing::Phased);
    /// assert_eq!(allele.position(), Some(0));
    /// ```
    pub fn position(&self) -> Option<usize> {
        self.position
    }

    /// Returns a mutable reference to the position of the allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::series::value::genotype::Phasing,
    ///     record_buf::samples::sample::value::genotype::Allele,
    /// };
    ///
    /// let mut allele = Allele::new(Some(0), Phasing::Phased);
    /// *allele.position_mut() = Some(1);
    /// assert_eq!(allele.position(), Some(1));
    /// ```
    pub fn position_mut(&mut self) -> &mut Option<usize> {
        &mut self.position
    }

    /// Returns the phasing of the allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::series::value::genotype::Phasing,
    ///     record_buf::samples::sample::value::genotype::Allele,
    /// };
    ///
    /// let allele = Allele::new(Some(0), Phasing::Phased);
    /// assert_eq!(allele.phasing(), Phasing::Phased);
    /// ```
    pub fn phasing(&self) -> Phasing {
        self.phasing
    }

    /// Returns a mutable reference to the phasing of the allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::variant::{
    ///     record::samples::series::value::genotype::Phasing,
    ///     record_buf::samples::sample::value::genotype::Allele,
    /// };
    ///
    /// let mut allele = Allele::new(Some(0), Phasing::Phased);
    /// *allele.phasing_mut() = Phasing::Unphased;
    /// assert_eq!(allele.phasing(), Phasing::Unphased);
    /// ```
    pub fn phasing_mut(&mut self) -> &mut Phasing {
        &mut self.phasing
    }
}

/// An error returned when a raw VCF record genotype value allele fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The position is invalid.
    InvalidPosition(num::ParseIntError),
    /// The phasing is invalid.
    InvalidPhasing,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidPosition(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidPosition(_) => f.write_str("invalid position"),
            Self::InvalidPhasing => f.write_str("invalid phasing"),
        }
    }
}

impl FromStr for Allele {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let phasing = parse_phasing(&s[..1])?;
        let position = parse_position(&s[1..])?;

        Ok(Allele::new(position, phasing))
    }
}

pub(super) fn parse_phasing(s: &str) -> Result<Phasing, ParseError> {
    const PHASED: &str = "|";
    const UNPHASED: &str = "/";

    match s {
        PHASED => Ok(Phasing::Phased),
        UNPHASED => Ok(Phasing::Unphased),
        _ => Err(ParseError::InvalidPhasing),
    }
}

pub(super) fn parse_position(s: &str) -> Result<Option<usize>, ParseError> {
    const MISSING: &str = ".";

    match s {
        MISSING => Ok(None),
        _ => s.parse().map(Some).map_err(ParseError::InvalidPosition),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("/.".parse(), Ok(Allele::new(None, Phasing::Unphased)));
        assert_eq!("/0".parse(), Ok(Allele::new(Some(0), Phasing::Unphased)));
        assert_eq!("/13".parse(), Ok(Allele::new(Some(13), Phasing::Unphased)));

        assert_eq!("".parse::<Allele>(), Err(ParseError::Empty));
        assert!(matches!(
            "/ndls".parse::<Allele>(),
            Err(ParseError::InvalidPosition(_))
        ));
        assert_eq!(":0".parse::<Allele>(), Err(ParseError::InvalidPhasing));
    }
}
