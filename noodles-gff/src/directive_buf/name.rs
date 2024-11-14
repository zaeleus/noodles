//! GFF directive name.

pub mod other;
mod standard;

pub use self::other::Other;
use self::standard::Standard;

pub(super) const GFF_VERSION: Name = Name::Standard(Standard::GffVersion);
pub(super) const SEQUENCE_REGION: Name = Name::Standard(Standard::SequenceRegion);
pub(super) const FEATURE_ONTOLOGY: Name = Name::Standard(Standard::FeatureOntology);
pub(super) const ATTRIBUTE_ONTOLOGY: Name = Name::Standard(Standard::AttributeOntology);
pub(super) const SOURCE_ONTOLOGY: Name = Name::Standard(Standard::SourceOntology);
pub(super) const SPECIES: Name = Name::Standard(Standard::Species);
pub(super) const GENOME_BUILD: Name = Name::Standard(Standard::GenomeBuild);
pub(super) const FORWARD_REFERENCES_ARE_RESOLVED: Name =
    Name::Standard(Standard::ForwardReferencesAreResolved);
pub(super) const START_OF_FASTA: Name = Name::Standard(Standard::StartOfFasta);

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) enum Name {
    Standard(Standard),
    Other(Other),
}

impl From<&str> for Name {
    fn from(s: &str) -> Self {
        match Standard::new(s) {
            Some(name) => Self::Standard(name),
            None => Self::Other(Other(s.into())),
        }
    }
}
