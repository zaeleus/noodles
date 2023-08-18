use std::{error, fmt, num};

use super::field::{parse_key, parse_value};
use crate::header::record::value::{
    map::{
        self,
        filter::{tag, Tag},
        Filter, OtherFields,
    },
    Map,
};

#[derive(Clone, Debug, Eq, PartialEq)]
enum ParseErrorKind {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    InvalidKey(super::field::key::ParseError),
    InvalidValue(super::field::value::ParseError),
    MissingId,
    MissingDescription,
    InvalidDescription,
    InvalidIdx(num::ParseIntError),
    InvalidOther(
        map::tag::Other<tag::Standard>,
        super::field::value::ParseError,
    ),
    DuplicateTag(Tag),
}

/// An error returned when a VCF header record filter map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError {
    id: Option<String>,
    kind: ParseErrorKind,
}

impl ParseError {
    fn new(id: Option<String>, kind: ParseErrorKind) -> Self {
        Self { id, kind }
    }

    pub(crate) fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match &self.kind {
            ParseErrorKind::InvalidMap(e) => Some(e),
            ParseErrorKind::InvalidField(e) => Some(e),
            ParseErrorKind::InvalidKey(e) => Some(e),
            ParseErrorKind::InvalidValue(e) => Some(e),
            ParseErrorKind::InvalidIdx(e) => Some(e),
            ParseErrorKind::InvalidOther(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ParseErrorKind::InvalidMap(_) => write!(f, "invalid map"),
            ParseErrorKind::InvalidField(_) => write!(f, "invalid field"),
            ParseErrorKind::InvalidKey(_) => write!(f, "invalid key"),
            ParseErrorKind::InvalidValue(_) => write!(f, "invalid value"),
            ParseErrorKind::MissingId => write!(f, "missing ID"),
            ParseErrorKind::MissingDescription => write!(f, "missing description"),
            ParseErrorKind::InvalidDescription => write!(f, "invalid description"),
            ParseErrorKind::InvalidIdx(_) => write!(f, "invalid IDX"),
            ParseErrorKind::InvalidOther(tag, _) => write!(f, "invalid other: {tag}"),
            ParseErrorKind::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_filter(src: &mut &[u8]) -> Result<(String, Map<Filter>), ParseError> {
    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut description = None;
    let mut idx = None;

    let mut other_fields = OtherFields::new();

    loop {
        let tag = parse_key(src)
            .map(Tag::from)
            .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidKey(e)))?;

        match tag {
            tag::ID => parse_id(src, &id).and_then(|v| try_replace(&mut id, &None, tag::ID, v))?,
            tag::DESCRIPTION => parse_description(src, &id)
                .and_then(|v| try_replace(&mut description, &id, tag::DESCRIPTION, v))?,
            tag::IDX => {
                parse_idx(src, &id).and_then(|v| try_replace(&mut idx, &id, tag::IDX, v))?
            }
            Tag::Other(t) => parse_other(src, &id, &t)
                .and_then(|value| try_insert(&mut other_fields, &id, t, value))?,
        }

        let has_separator = super::field::consume_separator(src)
            .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidField(e)))?;

        if !has_separator {
            break;
        }
    }

    super::consume_suffix(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidMap(e)))?;

    let id = id.ok_or_else(|| ParseError::new(None, ParseErrorKind::MissingId))?;
    let description = description
        .ok_or_else(|| ParseError::new(Some(id.clone()), ParseErrorKind::MissingDescription))?;

    Ok((
        id,
        Map {
            inner: Filter { description, idx },
            other_fields,
        },
    ))
}

fn parse_id(src: &mut &[u8], id: &Option<String>) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
}

fn parse_description(src: &mut &[u8], id: &Option<String>) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
}

fn parse_idx(src: &mut &[u8], id: &Option<String>) -> Result<usize, ParseError> {
    parse_value(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidIdx(e)))
        })
}

fn parse_other(
    src: &mut &[u8],
    id: &Option<String>,
    tag: &map::tag::Other<tag::Standard>,
) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidOther(tag.clone(), e)))
}

fn try_replace<T>(
    option: &mut Option<T>,
    id: &Option<String>,
    tag: Tag,
    value: T,
) -> Result<(), ParseError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(ParseError::new(
            id.clone(),
            ParseErrorKind::DuplicateTag(tag),
        ))
    }
}

fn try_insert(
    other_fields: &mut OtherFields<tag::Standard>,
    id: &Option<String>,
    tag: map::tag::Other<tag::Standard>,
    value: String,
) -> Result<(), ParseError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(entry) => {
            let (t, _) = entry.remove_entry();
            Err(ParseError::new(
                id.clone(),
                ParseErrorKind::DuplicateTag(Tag::Other(t)),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_filter() {
        let mut src = &br#"<ID=PASS,Description="All filters passed">"#[..];

        let id = String::from("PASS");
        let map = Map::<Filter>::pass();
        let expected = (id, map);

        assert_eq!(parse_filter(&mut src), Ok(expected));
    }
}
