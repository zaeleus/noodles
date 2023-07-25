use std::{error, fmt};

use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            program::{tag, Tag},
            OtherFields, Program,
        },
        Map,
    },
};

use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value};

/// An error returned when a SAM header program record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),

    MissingField(Tag),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(_) => write!(f, "invalid field"),
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_program(
    src: &mut &[u8],
    ctx: &Context,
) -> Result<(String, Map<Program>), ParseError> {
    let mut id = None;
    let mut name = None;
    let mut command_line = None;
    let mut previous_id = None;
    let mut description = None;
    let mut version = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidField)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::ID => parse_string(src).and_then(|v| try_replace(&mut id, ctx, tag::ID, v))?,
            tag::NAME => {
                parse_string(src).and_then(|v| try_replace(&mut name, ctx, tag::NAME, v))?
            }
            tag::COMMAND_LINE => parse_string(src)
                .and_then(|v| try_replace(&mut command_line, ctx, tag::COMMAND_LINE, v))?,
            tag::PREVIOUS_ID => parse_string(src)
                .and_then(|v| try_replace(&mut previous_id, ctx, tag::PREVIOUS_ID, v))?,
            tag::DESCRIPTION => parse_string(src)
                .and_then(|v| try_replace(&mut description, ctx, tag::DESCRIPTION, v))?,
            tag::VERSION => {
                parse_string(src).and_then(|v| try_replace(&mut version, ctx, tag::VERSION, v))?
            }
            Tag::Other(t) => parse_value(src)
                .map_err(ParseError::InvalidField)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let id = id.ok_or(ParseError::MissingField(tag::ID))?;

    Ok((
        id,
        Map {
            inner: Program {
                name,
                command_line,
                previous_id,
                description,
                version,
            },
            other_fields,
        },
    ))
}

fn parse_string(src: &mut &[u8]) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(ParseError::InvalidField)
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
    V: Into<String>,
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

        let expected = (String::from("pg0"), Map::<Program>::default());

        assert_eq!(actual, Ok(expected));
    }
}
