mod version;

use std::{error, fmt};

use bstr::{BStr, BString};

use self::version::parse_version;
use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value, value};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            header::{tag, Tag},
            tag::Other,
            Header, OtherFields,
        },
        Map,
    },
};

/// An error returned when a SAM header header record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),
    InvalidTag(super::field::tag::ParseError),
    InvalidValue(value::ParseError),
    MissingVersion,
    InvalidVersion(version::ParseError),
    InvalidOther(Other<tag::Standard>, value::ParseError),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidTag(e) => Some(e),
            Self::InvalidOther(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(_) => write!(f, "invalid field"),
            Self::InvalidTag(_) => write!(f, "invalid tag"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
            Self::MissingVersion => write!(f, "missing version ({}) field", tag::VERSION),
            Self::InvalidVersion(_) => write!(f, "invalid version ({})", tag::VERSION),
            Self::InvalidOther(tag, _) => write!(f, "invalid other ({tag})"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_header(src: &mut &[u8], ctx: &Context) -> Result<Map<Header>, ParseError> {
    let mut version = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::VERSION => {
                parse_value(src)
                    .map_err(ParseError::InvalidValue)
                    .and_then(|buf| parse_version(buf).map_err(ParseError::InvalidVersion))
                    .and_then(|v| try_replace(&mut version, ctx, tag::VERSION, v))?;
            }
            Tag::Other(t) => parse_other(src, t)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let version = version.ok_or(ParseError::MissingVersion)?;

    Ok(Map {
        inner: Header { version },
        other_fields,
    })
}

fn parse_other<'a>(src: &mut &'a [u8], tag: Other<tag::Standard>) -> Result<&'a BStr, ParseError> {
    parse_value(src).map_err(|e| ParseError::InvalidOther(tag, e))
}

fn try_replace<T>(
    option: &mut Option<T>,
    ctx: &Context,
    tag: Tag,
    value: T,
) -> Result<(), ParseError> {
    if option.replace(value).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(tag))
    } else {
        Ok(())
    }
}

fn try_insert<V>(
    other_fields: &mut OtherFields<tag::Standard>,
    ctx: &Context,
    tag: map::tag::Other<tag::Standard>,
    value: V,
) -> Result<(), ParseError>
where
    V: Into<BString>,
{
    if other_fields.insert(tag, value.into()).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(Tag::Other(tag)))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::record::value::map::header::Version;

    #[test]
    fn test_parse_header() {
        let mut src = &b"\tVN:1.6"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_header(&mut src, &ctx),
            Ok(Map::<Header>::new(Version::new(1, 6)))
        );
    }

    #[test]
    fn test_parse_header_with_missing_version() {
        let mut src = &b"\tSO:coordinate"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_header(&mut src, &ctx),
            Err(ParseError::MissingVersion)
        );
    }
}
