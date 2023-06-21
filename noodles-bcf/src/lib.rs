#![warn(missing_docs)]

//! **noodles-bcf** handles the reading and writing of the BCF format.

#[cfg(feature = "async")]
mod r#async;

pub mod header;
pub mod indexed_reader;
pub mod lazy;
pub mod reader;
pub(crate) mod record;
mod writer;

pub use self::{indexed_reader::IndexedReader, reader::Reader, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

static MAGIC_NUMBER: &[u8] = b"BCF";
