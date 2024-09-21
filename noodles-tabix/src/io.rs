//! Tabix I/O.

pub mod indexed_reader;
mod reader;
mod writer;

pub use self::{reader::Reader, writer::Writer};
