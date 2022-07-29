//! Async gzip index.
//!
//! A [gzip index] (GZI) is a list of compressed and uncompressed offset pairs for a gzipped file.
//!
//! [gzip index]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

mod reader;

pub use self::reader::Reader;
