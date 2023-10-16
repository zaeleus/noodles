use super::Record;

/// An immutable, lazily-evalulated GFF line.
pub enum Line {
    /// A directive (`##`).
    Directive(String),
    /// A comment (`#`),
    Comment(String),
    /// A record.
    Record(Record),
}

impl Default for Line {
    fn default() -> Self {
        Self::Comment(String::new())
    }
}

impl From<Line> for String {
    fn from(line: Line) -> Self {
        match line {
            Line::Directive(s) => s,
            Line::Comment(s) => s,
            Line::Record(record) => record.into(),
        }
    }
}
