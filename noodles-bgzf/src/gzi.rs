//! GZI index.

mod reader;

/// A GZI Index, representing pairs of compressed and uncompressed offsets in a BGZF file.
pub type Index = Vec<(u64, u64)>;