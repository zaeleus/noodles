mod kind;

use std::io;

pub use self::kind::Kind;
use super::Record;

const COMMENT_PREFIX: char = '#';
const COMMENT_START: usize = 1;

/// A GTF line.
#[derive(Clone, Eq, PartialEq)]
pub struct Line(pub(crate) String);

impl Line {
    /// Returns the kind of line.
    pub fn kind(&self) -> Kind {
        if self.0.starts_with(COMMENT_PREFIX) {
            Kind::Comment
        } else {
            Kind::Record
        }
    }

    /// Returns the line as a comment.
    pub fn as_comment(&self) -> Option<&str> {
        match self.kind() {
            Kind::Comment => Some(&self.0[COMMENT_START..]),
            Kind::Record => None,
        }
    }

    /// Returns the line as a record.
    pub fn as_record(&self) -> Option<io::Result<Record<'_>>> {
        match self.kind() {
            Kind::Comment => None,
            Kind::Record => Some(Record::try_new(&self.0)),
        }
    }
}

impl AsRef<str> for Line {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Default for Line {
    fn default() -> Self {
        Self(COMMENT_PREFIX.into())
    }
}
