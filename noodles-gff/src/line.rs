//! GFF line.

mod kind;

use std::io;

pub use self::kind::Kind;
use super::{Directive, Record};

const COMMENT_PREFIX: char = '#';

const DIRECTIVE_START: usize = 2;
const COMMENT_START: usize = 1;

/// A GFF line.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Line(pub(crate) String);

impl Line {
    /// Returns the kind of line.
    pub fn kind(&self) -> Kind {
        if let Some(src) = self.0.strip_prefix(COMMENT_PREFIX) {
            if src.starts_with(COMMENT_PREFIX) {
                Kind::Directive
            } else {
                Kind::Comment
            }
        } else {
            Kind::Record
        }
    }

    /// Returns the line as a directive.
    pub fn as_directive(&self) -> Option<Directive<'_>> {
        match self.kind() {
            Kind::Directive => Some(Directive::new(&self.0.as_bytes()[DIRECTIVE_START..])),
            _ => None,
        }
    }

    /// Returns the line as a comment.
    pub fn as_comment(&self) -> Option<&str> {
        match self.kind() {
            Kind::Comment => Some(&self.0[COMMENT_START..]),
            _ => None,
        }
    }

    /// Returns the line as a record.
    pub fn as_record(&self) -> Option<io::Result<Record<'_>>> {
        match self.kind() {
            Kind::Record => Some(Record::try_new(self.0.as_bytes())),
            _ => None,
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
        Self(String::from("#"))
    }
}
