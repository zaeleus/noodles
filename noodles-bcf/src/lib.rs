//! **noodles-bcf** handles the reading and writing of the BCF format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;

#[cfg(feature = "async")]
pub use self::r#async::io::{Reader as AsyncReader, Writer as AsyncWriter};

const MAGIC_NUMBER: [u8; 3] = *b"BCF";
