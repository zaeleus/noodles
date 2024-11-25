//! GFF directives.

pub mod genome_build;
pub mod gff_version;
pub mod key;
pub mod sequence_region;

pub use self::{
    genome_build::GenomeBuild, gff_version::GffVersion, sequence_region::SequenceRegion,
};

use std::{error, fmt, str::FromStr};

pub(crate) const PREFIX: &str = "##";

/// A GFF directive.
///
/// This is also called a pragma or metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DirectiveBuf {
    /// The GFF version (`gff-version`).
    GffVersion(GffVersion),
    /// A reference to a sequence segment (`sequence-region`).
    SequenceRegion(SequenceRegion),
    /// The genome build used for the start and end positions (`genome-build`).
    GenomeBuild(GenomeBuild),
    /// Any other directive.
    Other(String, Option<String>),
}

/// An error returned when a raw GFF directive fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The directive prefix (`##`) is missing.
    MissingPrefix,
    /// The directive name is missing.
    MissingKey,
    /// The directive value is missing.
    MissingValue,
    /// The GFF version is invalid.
    InvalidGffVersion(gff_version::ParseError),
    /// A sequence region is invalid.
    InvalidSequenceRegion(sequence_region::ParseError),
    /// A genome build is invalid.
    InvalidGenomeBuild(genome_build::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidGffVersion(e) => Some(e),
            Self::InvalidSequenceRegion(e) => Some(e),
            Self::InvalidGenomeBuild(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPrefix => f.write_str("directive prefix is missing"),
            Self::MissingKey => f.write_str("directive key is missing"),
            Self::MissingValue => f.write_str("directive value is missing"),
            Self::InvalidGffVersion(_) => f.write_str("invalid GFF version"),
            Self::InvalidSequenceRegion(_) => f.write_str("invalid sequence region"),
            Self::InvalidGenomeBuild(_) => f.write_str("invalid genome build"),
        }
    }
}

impl FromStr for DirectiveBuf {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !s.starts_with(PREFIX) {
            return Err(ParseError::MissingPrefix);
        }

        let mut components = s[PREFIX.len()..].splitn(2, |c: char| c.is_ascii_whitespace());

        let key = components.next().ok_or(ParseError::MissingKey)?;

        match key {
            key::GFF_VERSION => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidGffVersion))
                .map(Self::GffVersion),
            key::SEQUENCE_REGION => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidSequenceRegion))
                .map(Self::SequenceRegion),
            key::GENOME_BUILD => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidGenomeBuild))
                .map(Self::GenomeBuild),
            _ => {
                let value = components.next().map(String::from);
                Ok(Self::Other(key.into(), value))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(
            "##noodles".parse(),
            Ok(DirectiveBuf::Other(String::from("noodles"), None)),
        );

        assert_eq!(
            "##noodles gff".parse(),
            Ok(DirectiveBuf::Other(
                String::from("noodles"),
                Some(String::from("gff"))
            )),
        );
    }
}
