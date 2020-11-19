mod reader;
mod record;
mod writer;

pub use self::{reader::Reader, record::Record, writer::Writer};

/// A CRAM index.
pub type Index = Vec<Record>;
