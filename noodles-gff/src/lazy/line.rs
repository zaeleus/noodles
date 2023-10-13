/// An immutable, lazily-evalulated GFF line.
pub enum Line {
    /// A directive (`##`).
    Directive(String),
    /// A comment (`#`),
    Comment(String),
    /// A record.
    Record(String),
}

impl From<Line> for String {
    fn from(line: Line) -> Self {
        match line {
            Line::Directive(s) => s,
            Line::Comment(s) => s,
            Line::Record(s) => s,
        }
    }
}
