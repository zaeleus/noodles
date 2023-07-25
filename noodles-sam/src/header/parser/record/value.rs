mod comment;
mod map;

use std::{error, fmt};

use self::comment::parse_comment;
use crate::header::{parser::Context, record::Kind, Record};

/// An error returned when a SAM header record value fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The header is invalid.
    InvalidHeader(map::header::ParseError),
    /// The reference sequence is invalid.
    InvalidReferenceSequence(map::reference_sequence::ParseError),
    /// The comment is invalid.
    InvalidComment(comment::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidHeader(e) => Some(e),
            Self::InvalidReferenceSequence(e) => Some(e),
            Self::InvalidComment(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidHeader(_) => write!(f, "invalid header"),
            Self::InvalidReferenceSequence(_) => write!(f, "invalid reference sequence"),
            Self::InvalidComment(_) => write!(f, "invalid comment"),
        }
    }
}

pub(super) fn parse_value(
    src: &mut &[u8],
    ctx: &Context,
    kind: Kind,
) -> Result<Record, ParseError> {
    match kind {
        Kind::Header => map::parse_header(src, ctx)
            .map(Record::Header)
            .map_err(ParseError::InvalidHeader),
        Kind::ReferenceSequence => map::parse_reference_sequence(src, ctx)
            .map(|(name, map)| Record::ReferenceSequence(name, map))
            .map_err(ParseError::InvalidReferenceSequence),
        Kind::ReadGroup => todo!(),
        Kind::Program => todo!(),
        Kind::Comment => parse_comment(src)
            .map(Record::Comment)
            .map_err(ParseError::InvalidComment),
    }
}
