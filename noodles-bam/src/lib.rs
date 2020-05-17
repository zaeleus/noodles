pub mod bai;
pub mod reader;
pub mod record;
mod reference;
mod writer;

pub use self::{
    reader::{Reader, Records},
    record::Record,
    reference::Reference,
    writer::Writer,
};

pub static MAGIC_NUMBER: &[u8] = b"BAM\x01";
