#![warn(missing_docs)]

//! **noodles-cram** handles the reading and writing of the CRAM format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod codecs;
pub mod container;
pub mod crai;
pub mod file_definition;
pub mod fs;
mod huffman;
pub mod io;
mod num;
pub mod record;

pub use self::{container::Container, file_definition::FileDefinition, record::Record};

#[deprecated(since = "0.78.0", note = "Use `cram::container` instead.")]
pub use self::container as data_container;

#[deprecated(since = "0.78.0", note = "Use `cram::Container` instead.")]
pub use self::container::Container as DataContainer;

#[deprecated(since = "0.76.0", note = "Use `cram::fs::index` instead.")]
pub use self::fs::index;

#[cfg(feature = "async")]
#[deprecated(since = "0.69.0", note = "Use `cram::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.69.0", note = "Use `cram::r#async::io::Writer` instead.")]
pub use self::r#async::io::Writer as AsyncWriter;

const MAGIC_NUMBER: [u8; 4] = *b"CRAM";
