#![deny(missing_docs)]

//! **noodles-vcf** handles the reading and writing of the VCF format.

pub mod header;
mod reader;
pub mod record;
mod writer;

pub use self::{header::Header, reader::Reader, record::Record, writer::Writer};
