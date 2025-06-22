use std::{error, fmt};

use crate::header::record::value::map::{Tag, tag};

/// An error returned when a SAM header record field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
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

pub fn parse_tag<S>(src: &mut &[u8]) -> Result<Tag<S>, ParseError>
where
    S: tag::Standard,
{
    let Some((buf, rest)) = src.split_first_chunk() else {
        return Err(ParseError::UnexpectedEof);
    };

    let tag = Tag::<S>::from(*buf);

    *src = rest;

    Ok(tag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag() {
        use crate::header::record::value::map::header;

        let mut src = &b"VN"[..];
        assert_eq!(parse_tag(&mut src), Ok(header::tag::VERSION));

        let mut src = &b""[..];
        assert_eq!(
            parse_tag::<header::tag::Standard>(&mut src),
            Err(ParseError::UnexpectedEof)
        );

        let mut src = &b"V"[..];
        assert_eq!(
            parse_tag::<header::tag::Standard>(&mut src),
            Err(ParseError::UnexpectedEof)
        );
    }
}
