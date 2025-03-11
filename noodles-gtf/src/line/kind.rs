/// A GTF line kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    /// A comment (`#`).
    Comment,
    /// A record.
    Record,
}
