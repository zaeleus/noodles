//! GFF lines.

use bstr::BString;

use super::{feature::RecordBuf, DirectiveBuf};

/// A GFF line.
#[derive(Clone, Debug, PartialEq)]
pub enum LineBuf {
    /// A directive (`##`).
    Directive(DirectiveBuf),
    /// A comment (`#`),
    Comment(BString),
    /// A record.
    Record(RecordBuf),
}
