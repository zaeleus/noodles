//! VCF record genotype value allele.

pub mod phasing;

pub use self::phasing::Phasing;

use std::{error, fmt, num, str::FromStr};

const MISSING_POSITION: &str = ".";

/// A VCF record genotype value allele.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Allele {
    position: Option<usize>,
    phasing: Option<Phasing>,
}

impl Allele {
    /// Creates a VCF record genotype value allele.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::record::genotype::field::value::genotype::Allele;
    /// let allele = Allele::new(Some(0), None);
    /// ```
    pub fn new(position: Option<usize>, phasing: Option<Phasing>) -> Self {
        Self { position, phasing }
    }

    /// Returns the phasing of the allele.
    pub fn phasing(&self) -> Option<Phasing> {
        self.phasing
    }

    /// Returns a mutable reference to the phasing of the allele.
    pub fn phasing_mut(&mut self) -> &mut Option<Phasing> {
        &mut self.phasing
    }

    /// Returns the position of the allele.
    pub fn position(&self) -> Option<usize> {
        self.position
    }

    /// Returns a mutable reference to the position of the allele.
    pub fn position_mut(&mut self) -> &mut Option<usize> {
        &mut self.position
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
    InvalidPhasing(phasing::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidPosition(e) => write!(f, "invalid position: {}", e),
            Self::InvalidPhasing(e) => write!(f, "invalid phasing: {}", e),
        }
    }
}

impl FromStr for Allele {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        match s[..1].parse() {
            Ok(phasing) => {
                let position = parse_position(&s[1..])?;
                Ok(Allele::new(position, Some(phasing)))
            }
            Err(e) => {
                if let Ok(position) = parse_position(s) {
                    Ok(Allele::new(position, None))
                } else {
                    Err(ParseError::InvalidPhasing(e))
                }
            }
        }
    }
}

fn parse_position(s: &str) -> Result<Option<usize>, ParseError> {
    if s == MISSING_POSITION {
        Ok(None)
    } else {
        s.parse().map(Some).map_err(ParseError::InvalidPosition)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(".".parse(), Ok(Allele::new(None, None)));
        assert_eq!("0".parse(), Ok(Allele::new(Some(0), None)));
        assert_eq!("/.".parse(), Ok(Allele::new(None, Some(Phasing::Unphased))));
        assert_eq!(
            "/0".parse(),
            Ok(Allele::new(Some(0), Some(Phasing::Unphased)))
        );
        assert_eq!(
            "/13".parse(),
            Ok(Allele::new(Some(13), Some(Phasing::Unphased)))
        );

        assert_eq!("".parse::<Allele>(), Err(ParseError::Empty));
        assert!(matches!(
            "/ndls".parse::<Allele>(),
            Err(ParseError::InvalidPosition(_))
        ));
        assert!(matches!(
            ":0".parse::<Allele>(),
            Err(ParseError::InvalidPhasing(_))
        ));
    }
}
