//! Alignment format utilities.

#[cfg(feature = "async")]
pub mod r#async;

mod index;
pub mod io;
pub mod iter;
mod record;

pub use self::{index::Index, record::Record};
