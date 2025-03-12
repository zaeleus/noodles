use std::io;

use super::{Line, RecordBuf};

/// A GTF line.
#[derive(Clone, Debug, PartialEq)]
pub enum LineBuf {
    /// A comment (`#`).
    Comment(String),
    /// A record.
    Record(RecordBuf),
}

impl TryFrom<Line> for LineBuf {
    type Error = io::Error;

    fn try_from(line: Line) -> Result<Self, Self::Error> {
        if let Some(s) = line.as_comment() {
            Ok(Self::Comment(s.into()))
        } else if let Some(result) = line.as_record() {
            result.and_then(RecordBuf::try_from).map(Self::Record)
        } else {
            unreachable!()
        }
    }
}
