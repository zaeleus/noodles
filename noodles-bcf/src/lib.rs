//! **noodles-bcf** handles the reading and writing of the BCF format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;
