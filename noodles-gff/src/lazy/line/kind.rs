/// A GFF line kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    /// A directive (`##`).
    Directive,
    /// A comment (`#`).
    Comment,
    /// A record.
    Record,
}
