//! SAM header header group order.

/// Records are not grouped (`none`).
pub const NONE: &[u8] = b"none";

/// Records are grouped by name (`query`).
pub const QUERY: &[u8] = b"query";

/// Alignments are grouped by reference sequence and position (`reference`).
pub const REFERENCE: &[u8] = b"reference";
