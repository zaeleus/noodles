//! GTF line.

mod kind;

use std::io;

use bstr::{BStr, BString, ByteSlice};

pub use self::kind::Kind;
use super::Record;

const COMMENT_PREFIX: u8 = b'#';
const COMMENT_START: usize = 1;

/// A GTF line.
#[derive(Clone, Eq, PartialEq)]
pub struct Line(pub(crate) BString);

impl Line {
    /// Returns the kind of line.
    pub fn kind(&self) -> Kind {
        if self.0.starts_with(&[COMMENT_PREFIX]) {
            Kind::Comment
        } else {
            Kind::Record
        }
    }

    /// Returns the line as a comment.
    pub fn as_comment(&self) -> Option<&BStr> {
        match self.kind() {
            Kind::Comment => Some(self.0[COMMENT_START..].as_bstr()),
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

impl AsRef<BStr> for Line {
    fn as_ref(&self) -> &BStr {
        self.0.as_ref()
    }
}

impl Default for Line {
    fn default() -> Self {
        Self(BString::from("#"))
    }
}
