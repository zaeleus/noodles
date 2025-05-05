//! Alignment format utilities.

#[cfg(feature = "async")]
pub mod r#async;

pub mod io;
pub mod iter;
mod record;

pub use self::record::Record;
