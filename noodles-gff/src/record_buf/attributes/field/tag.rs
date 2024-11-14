//! GFF record attributes field tag.

/// The feature ID (`ID`).
pub const ID: &str = "ID";

/// The feature display name (`Name`).
pub const NAME: &str = "Name";

/// A secondary name for the feature (`Alias`).
pub const ALIAS: &str = "Alias";

/// The feature parent (`Parent`).
pub const PARENT: &str = "Parent";

/// The target alignment (`Tag`).
pub const TARGET: &str = "Target";

/// The alignment of the feature to the target (`Gap`).
pub const GAP: &str = "Gap";

/// Temporal relationship disambiguation (`Derives_from`).
pub const DERIVES_FROM: &str = "Derives_from";

/// A comment (`Note`).
pub const NOTE: &str = "Note";

/// A database cross-reference (`Dbxref`).
pub const DBXREF: &str = "Dbxref";

/// A cross-reference to an ontology term (`Ontology_term`).
pub const ONTOLOGY_TERM: &str = "Ontology_term";

/// Indicates whether a feature is circular (`Is_circular`).
pub const IS_CIRCULAR: &str = "Is_circular";

/// A GFF record attributes field tag.
pub type Tag = String;
