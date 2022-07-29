//! gzip index.
//!
//! A [gzip index] (GZI) is a list of compressed and uncompressed offset pairs for a gzipped file.
//!
//! [GZ index]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

mod reader;

pub use self::reader::Reader;

/// A gzip index.
pub type Index = Vec<(u64, u64)>;
