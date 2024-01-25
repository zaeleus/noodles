#![warn(missing_docs)]

//! **noodles-bcf** handles the reading and writing of the BCF format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod header;
pub mod io;
pub mod lazy;
pub(crate) mod record;

#[cfg(feature = "async")]
pub use self::r#async::io::Reader as AsyncReader;

static MAGIC_NUMBER: &[u8] = b"BCF";
