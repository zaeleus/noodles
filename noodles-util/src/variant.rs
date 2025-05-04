//! Variant format utilities.

#[cfg(feature = "async")]
pub mod r#async;

pub mod io;
mod record;

pub use self::record::Record;
