pub mod bai;
pub mod reader;
pub mod record;
mod writer;

pub use self::{
    reader::{Reader, Records},
    record::Record,
    writer::Writer,
};

pub static MAGIC_NUMBER: &[u8] = b"BAM\x01";
