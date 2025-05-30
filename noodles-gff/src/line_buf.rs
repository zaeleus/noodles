//! GFF lines.

use bstr::BString;

use super::{DirectiveBuf, feature::RecordBuf};

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
