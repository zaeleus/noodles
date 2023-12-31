//! Alignment record and fields.

mod any_record;
pub mod iter;
pub mod record;
pub mod record_buf;

pub use self::{any_record::AnyRecord, record_buf::RecordBuf};
