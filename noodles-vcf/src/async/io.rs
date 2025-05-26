//! Async VCF I/O.

mod reader;
mod writer;

#[deprecated(since = "0.79.0", note = "Use `vcf::r#async::io::Reader` instead.")]
pub use self::reader::Reader;

pub use self::writer::Writer;
