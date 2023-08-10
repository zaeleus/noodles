use std::{error, fmt, str};

/// An error returned when a VCF header record escaped string fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidEscapeSequence { b: u8 },
    InvalidUtf8(str::Utf8Error),
    UnexpectedEof,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidUtf8(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidEscapeSequence { b } => write!(f, "invalid escape sequence: '\\{b}'"),
            Self::InvalidUtf8(_) => write!(f, "invalid UTF-8"),
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub fn parse_raw_string<'a>(src: &mut &'a [u8]) -> Result<&'a str, ParseError> {
    const COMMA: u8 = b',';
    const GREATER_THAN_SIGN: u8 = b'>';

    if let Some(i) = src
        .iter()
        .position(|&b| matches!(b, COMMA | GREATER_THAN_SIGN))
    {
        let (buf, rest) = src.split_at(i);
        *src = rest;
        str::from_utf8(buf).map_err(ParseError::InvalidUtf8)
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

pub fn parse_escaped_string<'a>(src: &mut &'a [u8]) -> Result<&'a str, ParseError> {
    const BACKSLASH: u8 = b'\\';
    const QUOTATION_MARK: u8 = b'"';

    enum State {
        Normal,
        Escape,
        Done { offset: usize },
    }

    let mut state = State::Normal;

    for (i, &b) in src.iter().enumerate() {
        if b == BACKSLASH {
            state = State::Escape;
            continue;
        }

        match state {
            State::Normal => {
                if b == QUOTATION_MARK {
                    state = State::Done { offset: i };
                    break;
                }
            }
            State::Escape => {
                if matches!(b, BACKSLASH | QUOTATION_MARK) {
                    state = State::Normal;
                } else {
                    return Err(ParseError::InvalidEscapeSequence { b });
                }
            }
            State::Done { .. } => unreachable!(),
        }
    }

    if let State::Done { offset } = state {
        let (buf, rest) = src.split_at(offset);

        *src = &rest[1..];

        // FIXME: Unescaped escaped string.
        str::from_utf8(buf).map_err(ParseError::InvalidUtf8)
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_raw_string() {
        let mut src = &b"noodles-vcf,"[..];
        assert_eq!(parse_raw_string(&mut src), Ok("noodles-vcf"));

        let mut src = &b"noodles-vcf>"[..];
        assert_eq!(parse_raw_string(&mut src), Ok("noodles-vcf"));

        let mut src = &b"noodles-vcf"[..];
        assert_eq!(parse_raw_string(&mut src), Err(ParseError::UnexpectedEof));
    }

    #[test]
    fn test_parse_escaped_string() {
        let mut src = &br#"noodles-vcf","#[..];
        assert_eq!(parse_escaped_string(&mut src), Ok("noodles-vcf"));

        // FIXME
        // let mut src = &br#"noodles-\"vcf\"","#[..];
        // assert_eq!(parse_escaped_string(&mut src), Ok("noodles-\"vcf\""));

        // FIXME
        // let mut src = &br#"noodles\\vcf","#[..];
        // assert_eq!(parse_escaped_string(&mut src), Ok("noodles\\vcf"));

        let mut src = &br#"noodles\nvcf","#[..];
        assert_eq!(
            parse_escaped_string(&mut src),
            Err(ParseError::InvalidEscapeSequence { b: b'n' })
        );

        let mut src = &br#"noodles-vcf"#[..];
        assert_eq!(
            parse_escaped_string(&mut src),
            Err(ParseError::UnexpectedEof)
        );
    }
}
