mod comment;

use std::{error, fmt};

use self::comment::parse_comment;
use crate::header::{record::Kind, Record};

/// An error returned when a SAM header record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The comment is invalid.
    InvalidComment(comment::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidComment(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidComment(_) => write!(f, "invalid comment"),
        }
    }
}

pub(super) fn parse_value(src: &mut &[u8], kind: Kind) -> Result<Record, ParseError> {
    match kind {
        Kind::Header => todo!(),
        Kind::ReferenceSequence => todo!(),
        Kind::ReadGroup => todo!(),
        Kind::Program => todo!(),
        Kind::Comment => parse_comment(src)
            .map(Record::Comment)
            .map_err(ParseError::InvalidComment),
    }
}
