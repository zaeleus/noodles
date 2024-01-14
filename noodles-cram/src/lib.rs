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
pub use self::r#async::{Reader as AsyncReader, Writer as AsyncWriter};

static MAGIC_NUMBER: &[u8] = b"CRAM";
