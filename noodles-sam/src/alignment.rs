//! Alignment record and fields.

mod any_record;
pub mod iter;
pub mod record;

pub use self::record::Record;

#[doc(hidden)]
pub use self::any_record::AnyRecord;
