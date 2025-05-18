//! CRAM index.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;

/// A CRAM index.
pub type Index = Vec<Record>;
