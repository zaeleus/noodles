//! GFF lines.

use super::{feature::RecordBuf, DirectiveBuf};

/// A GFF line.
#[derive(Clone, Debug, PartialEq)]
pub enum LineBuf {
    /// A directive (`##`).
    Directive(DirectiveBuf),
    /// A comment (`#`),
    Comment(String),
    /// A record.
    Record(RecordBuf),
}
