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
    /// The ontology used for the feature types (`feature-ontology`).
    FeatureOntology(String),
    /// The ontology used for the attributes (`attribute-ontology`).
    AttributeOntology(String),
    /// The ontology used for the sources (`source-ontology`).
    SourceOntology(String),
    /// The species the annotations apply to (`species`).
    Species(String),
    /// The genome build used for the start and end positions (`genome-build`).
    GenomeBuild(GenomeBuild),
    /// A marker indicating that all forward references to feature IDs have been resolved (`#`).
    ForwardReferencesAreResolved,
    /// A marker indicating the end of the records list and start of a bundled reference sequences
    /// (`FASTA`).
    StartOfFasta,
    /// A nonstandard directive.
    Other(String, Option<String>),
}

impl fmt::Display for DirectiveBuf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GffVersion(version) => write!(f, "{PREFIX}{} {version}", key::GFF_VERSION),
            Self::SequenceRegion(sequence_region) => {
                write!(f, "{PREFIX}{} {sequence_region}", key::SEQUENCE_REGION)
            }
            Self::FeatureOntology(uri) => write!(f, "{PREFIX}{} {uri}", key::FEATURE_ONTOLOGY),
            Self::AttributeOntology(uri) => write!(f, "{PREFIX}{} {uri}", key::ATTRIBUTE_ONTOLOGY),
            Self::SourceOntology(uri) => write!(f, "{PREFIX}{} {uri}", key::SOURCE_ONTOLOGY),
            Self::Species(uri) => write!(f, "{PREFIX}{} {uri}", key::SPECIES),
            Self::GenomeBuild(genome_build) => {
                write!(f, "{PREFIX}{} {genome_build}", key::GENOME_BUILD)
            }
            Self::ForwardReferencesAreResolved => {
                write!(f, "{PREFIX}{}", key::FORWARD_REFERENCES_ARE_RESOLVED)
            }
            Self::StartOfFasta => write!(f, "{PREFIX}{}", key::START_OF_FASTA),
            Self::Other(key, value) => {
                write!(f, "{PREFIX}{key}")?;

                if let Some(v) = value {
                    write!(f, " {v}")?;
                }

                Ok(())
            }
        }
    }
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
            key::FEATURE_ONTOLOGY => components
                .next()
                .map(|s| Self::FeatureOntology(s.into()))
                .ok_or(ParseError::MissingValue),
            key::ATTRIBUTE_ONTOLOGY => components
                .next()
                .map(|s| Self::AttributeOntology(s.into()))
                .ok_or(ParseError::MissingValue),
            key::SOURCE_ONTOLOGY => components
                .next()
                .map(|s| Self::SourceOntology(s.into()))
                .ok_or(ParseError::MissingValue),
            key::SPECIES => components
                .next()
                .map(|s| Self::Species(s.into()))
                .ok_or(ParseError::MissingValue),
            key::GENOME_BUILD => components
                .next()
                .ok_or(ParseError::MissingValue)
                .and_then(|s| s.parse().map_err(ParseError::InvalidGenomeBuild))
                .map(Self::GenomeBuild),
            key::FORWARD_REFERENCES_ARE_RESOLVED => Ok(Self::ForwardReferencesAreResolved),
            key::START_OF_FASTA => Ok(Self::StartOfFasta),
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

    #[test]
    fn test_fmt() {
        assert_eq!(
            DirectiveBuf::GffVersion(GffVersion::default()).to_string(),
            "##gff-version 3"
        );

        let directive =
            DirectiveBuf::SequenceRegion(SequenceRegion::new(String::from("sq0"), 8, 13));
        assert_eq!(directive.to_string(), "##sequence-region sq0 8 13");

        assert_eq!(
            DirectiveBuf::FeatureOntology(String::from("https://example.com/fo.obo")).to_string(),
            "##feature-ontology https://example.com/fo.obo"
        );

        assert_eq!(
            DirectiveBuf::AttributeOntology(String::from("https://example.com/ao.obo")).to_string(),
            "##attribute-ontology https://example.com/ao.obo"
        );

        assert_eq!(
            DirectiveBuf::SourceOntology(String::from("https://example.com/so.obo")).to_string(),
            "##source-ontology https://example.com/so.obo"
        );

        assert_eq!(
            DirectiveBuf::Species(String::from("https://example.com/species?id=1")).to_string(),
            "##species https://example.com/species?id=1"
        );

        let directive =
            DirectiveBuf::GenomeBuild(GenomeBuild::new(String::from("NDLS"), String::from("r1")));
        assert_eq!(directive.to_string(), "##genome-build NDLS r1");

        assert_eq!(
            DirectiveBuf::ForwardReferencesAreResolved.to_string(),
            "###"
        );
        assert_eq!(DirectiveBuf::StartOfFasta.to_string(), "##FASTA");

        let directive = DirectiveBuf::Other(String::from("noodles"), None);
        assert_eq!(directive.to_string(), "##noodles");

        let directive = DirectiveBuf::Other(String::from("noodles"), Some(String::from("gff")));
        assert_eq!(directive.to_string(), "##noodles gff");
    }
}
