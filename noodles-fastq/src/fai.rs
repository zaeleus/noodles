//! FASTQ index.

pub mod io;
mod record;

pub use self::record::Record;

/// A FASTQ index.
pub type Index = Vec<Record>;
