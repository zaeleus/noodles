mod reader;
mod record;

pub use self::{reader::Reader, record::Record};

/// A CRAM index.
pub type Index = Vec<Record>;
