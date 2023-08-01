use std::{error, fmt};

use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value, value};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            header::{
                group_order, sort_order, subsort_order, tag, version, GroupOrder, SortOrder,
                SubsortOrder, Tag, Version,
            },
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
    InvalidSortOrder(sort_order::ParseError),
    InvalidGroupOrder(group_order::ParseError),
    InvalidSubsortOrder(subsort_order::ParseError),
    InvalidOther(Other<tag::Standard>, value::ParseError),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidTag(e) => Some(e),
            Self::InvalidVersion(e) => Some(e),
            Self::InvalidSortOrder(e) => Some(e),
            Self::InvalidGroupOrder(e) => Some(e),
            Self::InvalidSubsortOrder(e) => Some(e),
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
            Self::InvalidSortOrder(_) => write!(f, "invalid sort order ({})", tag::SORT_ORDER),
            Self::InvalidGroupOrder(_) => write!(f, "invalid group order ({})", tag::GROUP_ORDER),
            Self::InvalidSubsortOrder(_) => {
                write!(f, "invalid subsort order ({})", tag::SUBSORT_ORDER)
            }
            Self::InvalidOther(tag, _) => write!(f, "invalid other ({tag})"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_header(src: &mut &[u8], ctx: &Context) -> Result<Map<Header>, ParseError> {
    let mut version = None;
    let mut sort_order = None;
    let mut group_order = None;
    let mut subsort_order = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

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
            Tag::Other(t) => parse_other(src, t)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let version = version.ok_or(ParseError::MissingVersion)?;

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

fn parse_version(src: &mut &[u8]) -> Result<Version, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidVersion))
}

fn parse_sort_order(src: &mut &[u8]) -> Result<SortOrder, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidSortOrder))
}

fn parse_group_order(src: &mut &[u8]) -> Result<GroupOrder, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidGroupOrder))
}

fn parse_subsort_order(src: &mut &[u8]) -> Result<SubsortOrder, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidSubsortOrder))
}

fn parse_other(src: &mut &[u8], tag: Other<tag::Standard>) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::InvalidOther(tag, e))
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
