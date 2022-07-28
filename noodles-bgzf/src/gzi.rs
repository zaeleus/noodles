//! A module for [GZI]. A GZ index contains pairs of compressed and uncompressed offsets in a
//! BGZF file. Values in the index are stored as little-endian 64-bit unsigned integers.
//!
//! [GZI]: http://www.htslib.org/doc/bgzip.html#GZI_FORMAT

#[cfg(feature = "async")]
pub mod r#async;

mod reader;

pub use self::reader::Reader;

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

/// A GZ Index represents pairs of compressed and uncompressed offsets in a BGZF file.
pub type Index = Vec<(u64, u64)>;
