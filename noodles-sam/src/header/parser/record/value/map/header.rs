use std::{error, fmt, str};

use indexmap::IndexMap;

use super::{parse_string, parse_tag};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            header::{
                group_order, sort_order, subsort_order, tag, version, GroupOrder, SortOrder,
                SubsortOrder, Tag, Version,
            },
            Header, OtherFields,
        },
        Map,
    },
};

/// An error returned when a SAM header header record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidDelimiter,
    InvalidTag(super::tag::ParseError),
    InvalidSeparator,
    InvalidVersion(version::ParseError),
    InvalidSortOrder(sort_order::ParseError),
    InvalidGroupOrder(group_order::ParseError),
    InvalidSubsortOrder(subsort_order::ParseError),
    InvalidValue(str::Utf8Error),
    MissingField(Tag),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidTag(e) => Some(e),
            Self::InvalidVersion(e) => Some(e),
            Self::InvalidSortOrder(e) => Some(e),
            Self::InvalidGroupOrder(e) => Some(e),
            Self::InvalidSubsortOrder(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDelimiter => write!(f, "invalid delimiter"),
            Self::InvalidTag(_) => write!(f, "invalid tag"),
            Self::InvalidSeparator => write!(f, "invalid separator"),
            Self::InvalidVersion(_) => write!(f, "invalid version"),
            Self::InvalidSortOrder(_) => write!(f, "invalid sort order"),
            Self::InvalidGroupOrder(_) => write!(f, "invalid group order"),
            Self::InvalidSubsortOrder(_) => write!(f, "invalid subsort order"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_header(src: &mut &[u8], ctx: &Context) -> Result<Map<Header>, ParseError> {
    let mut version = None;
    let mut sort_order = None;
    let mut group_order = None;
    let mut subsort_order = None;

    let mut other_fields = IndexMap::new();

    while !src.is_empty() {
        consume_delimiter(src)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src)?;

        match tag {
            tag::VERSION => {
                parse_version(src).and_then(|v| try_replace(&mut version, ctx, tag::VERSION, v))?
            }
            tag::SORT_ORDER => parse_sort_order(src)
                .and_then(|v| try_replace(&mut sort_order, ctx, tag::SORT_ORDER, v))?,
            tag::GROUP_ORDER => parse_group_order(src)
                .and_then(|v| try_replace(&mut group_order, ctx, tag::GROUP_ORDER, v))?,
            tag::SUBSORT_ORDER => parse_subsort_order(src)
                .and_then(|v| try_replace(&mut subsort_order, ctx, tag::SUBSORT_ORDER, v))?,
            Tag::Other(t) => parse_string(src)
                .map_err(ParseError::InvalidValue)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let version = version.ok_or(ParseError::MissingField(tag::VERSION))?;

    Ok(Map {
        inner: Header {
            version,
            sort_order,
            group_order,
            subsort_order,
        },
        other_fields,
    })
}

fn consume_delimiter(src: &mut &[u8]) -> Result<(), ParseError> {
    const PREFIX: u8 = b'\t';

    if let Some((b, rest)) = src.split_first() {
        if *b == PREFIX {
            *src = rest;
            return Ok(());
        }
    }

    Err(ParseError::InvalidDelimiter)
}

fn consume_separator(src: &mut &[u8]) -> Result<(), ParseError> {
    const SEPARATOR: u8 = b':';

    if let Some((b, rest)) = src.split_first() {
        if *b == SEPARATOR {
            *src = rest;
            return Ok(());
        }
    }

    Err(ParseError::InvalidSeparator)
}

fn parse_version(src: &mut &[u8]) -> Result<Version, ParseError> {
    parse_string(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidVersion))
}

fn parse_sort_order(src: &mut &[u8]) -> Result<SortOrder, ParseError> {
    parse_string(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidSortOrder))
}

fn parse_group_order(src: &mut &[u8]) -> Result<GroupOrder, ParseError> {
    parse_string(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidGroupOrder))
}

fn parse_subsort_order(src: &mut &[u8]) -> Result<SubsortOrder, ParseError> {
    parse_string(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidSubsortOrder))
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
    fn test_parse_header() -> Result<(), Box<dyn std::error::Error>> {
        let mut src = &b"\tVN:1.6\tSO:coordinate"[..];

        let ctx = Context::default();

        let actual = parse_header(&mut src, &ctx)?;
        let expected = Map::<Header>::builder()
            .set_version(Version::new(1, 6))
            .set_sort_order(SortOrder::Coordinate)
            .build()?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
