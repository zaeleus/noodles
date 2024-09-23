//! FASTQ index.

pub mod io;
mod record;

#[deprecated(since = "0.15.0", note = "Use `fai::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.15.0", note = "Use `fai::io::Writer` instead.")]
pub use self::io::Writer;

pub use self::record::Record;

/// A FASTQ index.
pub type Index = Vec<Record>;
