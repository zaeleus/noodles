//! Sequence repository adapters.

mod empty;
mod indexed_reader;
mod records;

pub use self::{empty::Empty, indexed_reader::IndexedReader};
