//! Alignment record and fields.

pub mod io;
pub mod iter;
pub mod record;
pub mod record_buf;

pub use self::{record::Record, record_buf::RecordBuf};
