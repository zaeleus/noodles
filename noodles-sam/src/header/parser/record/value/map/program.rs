use std::{error, fmt};

use bstr::{BStr, BString};

use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value, value};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            program::{tag, Tag},
            tag::Other,
            OtherFields, Program,
        },
        Map,
    },
};

/// An error returned when a SAM header program record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),
    InvalidTag(super::field::tag::ParseError),
    InvalidValue(value::ParseError),
    MissingId,
    InvalidId(value::ParseError),
    InvalidOther(Other<tag::Standard>, value::ParseError),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidTag(e) => Some(e),
            Self::InvalidId(e) => Some(e),
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
            Self::MissingId => write!(f, "missing ID field"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::InvalidOther(tag, _) => write!(f, "invalid other ({tag})"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_program(
    src: &mut &[u8],
    ctx: &Context,
) -> Result<(Vec<u8>, Map<Program>), ParseError> {
    let mut id = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::ID => parse_id(src).and_then(|v| try_replace(&mut id, ctx, tag::ID, v))?,
            Tag::Other(t) => parse_other(src, t)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let id = id.ok_or(ParseError::MissingId)?;

    Ok((
        id.to_vec(),
        Map {
            inner: Program,
            other_fields,
        },
    ))
}

fn parse_id<'a>(src: &mut &'a [u8]) -> Result<&'a BStr, ParseError> {
    parse_value(src).map_err(ParseError::InvalidId)
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

    #[test]
    fn test_parse_program() {
        let mut src = &b"\tID:pg0"[..];
        let ctx = Context::default();
        let actual = parse_program(&mut src, &ctx);
        let expected = (Vec::from("pg0"), Map::<Program>::default());
        assert_eq!(actual, Ok(expected));
    }

    #[test]
    fn test_parse_program_with_missing_id() {
        let mut src = &b"\tPN:pg0"[..];
        let ctx = Context::default();
        assert_eq!(parse_program(&mut src, &ctx), Err(ParseError::MissingId));
    }
}
