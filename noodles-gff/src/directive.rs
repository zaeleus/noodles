//! GFF directives.

use std::str::FromStr;

pub(crate) const PREFIX: &str = "##";

/// A GFF directive.
///
/// This is also called a pragma or metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Directive {
    /// The GFF version (`gff-version`).
    GffVersion(String),
    /// A reference to a sequence segment (`sequence-region`).
    SequenceRegion(String),
    /// The ontology used for the feature types (`feature-ontology`).
    FeatureOntology(String),
    /// The ontology used for the attributes (`attribute-ontology`).
    AttributeOntology(String),
    /// The ontology used for the sources (`source-ontology`).
    SourceOntology(String),
    /// The species the annotations apply to (`species`).
    Species(String),
    /// The genome build used for the start and end positions (`genome-build`).
    GenomeBuild(String),
    /// A marker indicating that all forward references to feature IDs have been resolved (`#`).
    ForwardReferencesAreResolved,
    /// A marker indicating the end of the records list and start of a bundled reference sequences
    /// (`FASTA`).
    StartOfFasta,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    Invalid(String),
    MissingPrefix,
    MissingName,
    MissingValue,
}

impl FromStr for Directive {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !s.starts_with(PREFIX) {
            return Err(ParseError::MissingPrefix);
        }

        let mut components = s[PREFIX.len()..].splitn(2, |c: char| c.is_ascii_whitespace());

        let name = components.next().ok_or_else(|| ParseError::MissingName)?;

        return match name {
            "gff-version" => components
                .next()
                .map(|s| Self::GffVersion(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "sequence-region" => components
                .next()
                .map(|s| Self::SequenceRegion(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "feature-ontology" => components
                .next()
                .map(|s| Self::FeatureOntology(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "attribute-ontology" => components
                .next()
                .map(|s| Self::AttributeOntology(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "source-ontology" => components
                .next()
                .map(|s| Self::SourceOntology(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "species" => components
                .next()
                .map(|s| Self::Species(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "genome-build" => components
                .next()
                .map(|s| Self::GenomeBuild(s.into()))
                .ok_or_else(|| ParseError::MissingValue),
            "#" => Ok(Self::ForwardReferencesAreResolved),
            "FASTA" => Ok(Self::StartOfFasta),
            _ => Err(ParseError::Invalid(name.into())),
        };
    }
}
