//! GFF directive key.

/// The GFF version (`gff-version`).
pub const GFF_VERSION: &str = "gff-version";

/// A reference to a sequence segment (`sequence-region`).
pub const SEQUENCE_REGION: &str = "sequence-region";

/// The ontology used for the feature types (`feature-ontology`).
pub const FEATURE_ONTOLOGY: &str = "feature-ontology";

/// The ontology used for the attributes (`attribute-ontology`).
pub const ATTRIBUTE_ONTOLOGY: &str = "attribute-ontology";

/// The ontology used for the sources (`source-ontology`).
pub const SOURCE_ONTOLOGY: &str = "source-ontology";

/// The species the annotations apply to (`species`).
pub const SPECIES: &str = "species";

/// The genome build used for the start and end positions (`genome-build`).
pub const GENOME_BUILD: &str = "genome-build";

/// A marker indicating that all forward references to feature IDs have been resolved (`#`).
pub const FORWARD_REFERENCES_ARE_RESOLVED: &str = "#";

/// A marker indicating the end of the records list and start of a bundled reference sequences
/// (`FASTA`).
pub const START_OF_FASTA: &str = "FASTA";
