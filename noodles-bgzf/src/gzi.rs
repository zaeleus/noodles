//! GZI index.

mod reader;

/// A GZI Index containing a number of entries, representing pairs of compressed and uncompressed offsets in a BGZF file.
pub struct Index {
  number_entries: u64,
  offsets: Vec<(u64, u64)>
}