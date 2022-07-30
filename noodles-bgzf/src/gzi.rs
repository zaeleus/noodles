//! gzip index.
//!
//! A [gzip index] (GZI) is a list of compressed and uncompressed offset pairs for a gzipped file.
//!
//! [gzip index]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

#[cfg(feature = "async")]
mod r#async;

mod reader;

pub use self::reader::Reader;

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

/// A gzip index.
pub type Index = Vec<(u64, u64)>;
