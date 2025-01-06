#![warn(missing_docs)]

//! **noodles-cram** handles the reading and writing of the CRAM format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod codecs;
pub(crate) mod container;
pub mod crai;
pub mod data_container;
pub mod file_definition;
mod huffman;
mod indexer;
pub mod io;
mod num;
pub mod record;

pub use self::{
    data_container::DataContainer, file_definition::FileDefinition, indexer::index, record::Record,
};

#[cfg(feature = "async")]
#[deprecated(since = "0.69.0", note = "Use `cram::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.69.0", note = "Use `cram::r#async::io::Writer` instead.")]
pub use self::r#async::io::Writer as AsyncWriter;

const MAGIC_NUMBER: [u8; 4] = *b"CRAM";
