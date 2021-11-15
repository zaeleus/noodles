#![warn(missing_docs)]

//! **noodles-bcf** handles the reading and writing of the BCF format.

#[cfg(feature = "async")]
mod r#async;

pub mod header;
pub mod reader;
pub mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

static MAGIC_NUMBER: &[u8] = b"BCF";
