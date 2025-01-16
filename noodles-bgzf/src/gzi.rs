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

#[deprecated(since = "0.35.0", note = "Use `gzi::fs::read` instead.")]
pub use self::fs::read;

#[deprecated(since = "0.35.0", note = "Use `gzi::io::Reader` instead.")]
pub use self::io::Reader;

#[cfg(feature = "async")]
#[deprecated(since = "0.35.0", note = "Use `bgzf::gzi::r#async::Reader` instead.")]
pub use self::r#async::Reader as AsyncReader;
