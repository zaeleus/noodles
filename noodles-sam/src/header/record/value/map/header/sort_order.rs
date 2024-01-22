//! SAM header header sort order.

/// The record order is unknown (`unknown`).
pub const UNKNOWN: &[u8] = b"unknown";

/// Records are not sorted (`unsorted`).
pub const UNSORTED: &[u8] = b"unsorted";

/// Records are sorted by name (`queryname`).
pub const QUERY_NAME: &[u8] = b"queryname";

/// Records are sorted by reference sequence and position (`coordinate`).
pub const COORDINATE: &[u8] = b"coordinate";
