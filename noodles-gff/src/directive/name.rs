mod other;
mod standard;

pub const GFF_VERSION: Name = Name::Standard(Standard::GffVersion);
pub const SEQUENCE_REGION: Name = Name::Standard(Standard::SequenceRegion);
pub const FEATURE_ONTOLOGY: Name = Name::Standard(Standard::FeatureOntology);
pub const ATTRIBUTE_ONTOLOGY: Name = Name::Standard(Standard::AttributeOntology);
pub const SOURCE_ONTOLOGY: Name = Name::Standard(Standard::SourceOntology);
pub const SPECIES: Name = Name::Standard(Standard::Species);
pub const GENOME_BUILD: Name = Name::Standard(Standard::GenomeBuild);
pub const FORWARD_REFERENCES_ARE_RESOLVED: Name =
    Name::Standard(Standard::ForwardReferencesAreResolved);
pub const START_OF_FASTA: Name = Name::Standard(Standard::StartOfFasta);

use self::{other::Other, standard::Standard};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Name {
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
