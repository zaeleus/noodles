use std::fmt;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Standard {
    GffVersion,
    SequenceRegion,
    FeatureOntology,
    AttributeOntology,
    SourceOntology,
    Species,
    GenomeBuild,
    ForwardReferencesAreResolved,
    StartOfFasta,
}

impl Standard {
    pub fn new(s: &str) -> Option<Self> {
        match s {
            "gff-version" => Some(Self::GffVersion),
            "sequence-region" => Some(Self::SequenceRegion),
            "feature-ontology" => Some(Self::FeatureOntology),
            "attribute-ontology" => Some(Self::AttributeOntology),
            "source-ontology" => Some(Self::SourceOntology),
            "species" => Some(Self::Species),
            "genome-build" => Some(Self::GenomeBuild),
            "#" => Some(Self::ForwardReferencesAreResolved),
            "FASTA" => Some(Self::StartOfFasta),
            _ => None,
        }
    }
}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::GffVersion => "gff-version",
            Self::SequenceRegion => "sequence-region",
            Self::FeatureOntology => "feature-ontology",
            Self::AttributeOntology => "attribute-ontology",
            Self::SourceOntology => "source-ontology",
            Self::Species => "species",
            Self::GenomeBuild => "genome-build",
            Self::ForwardReferencesAreResolved => "#",
            Self::StartOfFasta => "FASTA",
        }
    }
}

impl fmt::Display for Standard {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}
