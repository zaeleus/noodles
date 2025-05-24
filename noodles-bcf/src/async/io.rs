//! Async BCF I/O.

mod reader;
mod writer;

#[cfg(feature = "async")]
#[deprecated(since = "0.76.0", note = "Use `bcf::r#async::io::Reader` instead.")]
pub use self::reader::Reader;

#[deprecated(since = "0.76.0", note = "Use `bcf::r#async::io::Writer` instead.")]
pub use self::writer::Writer;
