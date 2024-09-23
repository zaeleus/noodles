//! FASTQ index.

pub mod io;
mod record;
mod writer;

#[deprecated(since = "0.15.0", note = "Use `fai::io::Reader` instead.")]
pub use self::io::Reader;

pub use self::{record::Record, writer::Writer};

/// A FASTQ index.
pub type Index = Vec<Record>;
