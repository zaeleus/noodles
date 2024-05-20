//! Async FASTA.

pub mod io;

#[deprecated(since = "0.39.0", note = "Use `gff::r#async::io::Reader` instead.")]
pub use self::io::Reader;
