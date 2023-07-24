use std::{error, fmt};

use crate::header::record::value::map::{tag, Tag};

/// An error returned when a SAM header record map value field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub(super) fn parse_tag<S>(src: &mut &[u8]) -> Result<Tag<S>, ParseError>
where
    S: tag::Standard,
{
    const TAG_LENGTH: usize = 2;

    if src.len() < TAG_LENGTH {
        return Err(ParseError::UnexpectedEof);
    }

    let (raw_tag, rest) = src.split_at(TAG_LENGTH);

    // SAFETY: `raw_tag` is `TAG_LENGTH`.
    let buf: [u8; TAG_LENGTH] = raw_tag.try_into().unwrap();
    let tag = Tag::<S>::from(buf);

    *src = rest;

    Ok(tag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() {
        use crate::header::record::value::map::header::tag;

        let mut src = &b"VN"[..];
        assert_eq!(parse_tag(&mut src), Ok(tag::VERSION));

        let mut src = &b""[..];
        assert_eq!(
            parse_tag::<tag::Standard>(&mut src),
            Err(ParseError::UnexpectedEof)
        );

        let mut src = &b"V"[..];
        assert_eq!(
            parse_tag::<tag::Standard>(&mut src),
            Err(ParseError::UnexpectedEof)
        );
    }
}
