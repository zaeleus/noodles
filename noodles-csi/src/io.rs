//! CSI I/O.

mod filter_by_region;
mod indexed_reader;
mod indexed_record;
mod indexed_records;
mod query;
pub mod reader;
pub(crate) mod writer;

pub use self::{
    filter_by_region::FilterByRegion, indexed_reader::IndexedReader, indexed_record::IndexedRecord,
    indexed_records::IndexedRecords, query::Query, reader::Reader, writer::Writer,
};

pub(crate) const MAGIC_NUMBER: [u8; 4] = *b"CSI\x01";
