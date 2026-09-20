//! Variant format utilities.

#[cfg(feature = "async")]
pub mod r#async;

pub mod fs;
mod index;
pub mod io;
mod record;

pub use self::{index::Index, record::Record};
