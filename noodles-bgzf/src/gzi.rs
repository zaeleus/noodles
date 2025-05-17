//! gzip index.
//!
//! A [gzip index] (GZI) is a list of compressed and uncompressed offset pairs for a gzipped file.
//!
//! [gzip index]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
mod index;
pub mod io;

pub use self::index::Index;
