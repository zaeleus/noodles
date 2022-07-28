//! A module for [GZI]. A GZ index contains pairs of compressed and uncompressed offsets in a
//! BGZF file. Values in the index are stored as little-endian 64-bit unsigned integers.
//!
//! [GZI]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

mod reader;

pub use self::reader::Reader;

/// A GZ Index represents pairs of compressed and uncompressed offsets in a BGZF file.
pub type Index = Vec<(u64, u64)>;
