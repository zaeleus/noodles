//! Sequence repository adapters.

mod empty;
mod indexed_reader;

pub use self::{empty::Empty, indexed_reader::IndexedReader};
