//! GFF directive key.

/// The GFF version (`gff-version`).
pub const GFF_VERSION: &[u8] = b"gff-version";

/// A reference to a sequence segment (`sequence-region`).
pub const SEQUENCE_REGION: &[u8] = b"sequence-region";

/// The ontology used for the feature types (`feature-ontology`).
pub const FEATURE_ONTOLOGY: &[u8] = b"feature-ontology";

/// The ontology used for the attributes (`attribute-ontology`).
pub const ATTRIBUTE_ONTOLOGY: &[u8] = b"attribute-ontology";

/// The ontology used for the sources (`source-ontology`).
pub const SOURCE_ONTOLOGY: &[u8] = b"source-ontology";

/// The species the annotations apply to (`species`).
pub const SPECIES: &[u8] = b"species";

/// The genome build used for the start and end positions (`genome-build`).
pub const GENOME_BUILD: &[u8] = b"genome-build";

/// A marker indicating that all forward references to feature IDs have been resolved (`#`).
pub const FORWARD_REFERENCES_ARE_RESOLVED: &[u8] = b"#";

/// A marker indicating the end of the records list and start of a bundled reference sequences
/// (`FASTA`).
pub const FASTA: &[u8] = b"FASTA";
