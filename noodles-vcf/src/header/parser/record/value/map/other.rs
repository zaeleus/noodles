use std::{error, fmt, str};

use super::field::{parse_key, parse_value};
use crate::header::record::value::{
    map::{
        self,
        other::{tag, Tag},
        Other, OtherFields,
    },
    Map,
};

/// An error returned when a VCF header record other map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    InvalidKey(super::field::key::ParseError),
    InvalidValue(super::field::value::ParseError),
    MissingId,
    InvalidValues(str::Utf8Error),
    InvalidOther(
        map::tag::Other<tag::StandardTag>,
        super::field::value::ParseError,
    ),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            ParseError::InvalidMap(e) => Some(e),
            ParseError::InvalidField(e) => Some(e),
            ParseError::InvalidKey(e) => Some(e),
            ParseError::InvalidValue(e) => Some(e),
            ParseError::InvalidValues(e) => Some(e),
            ParseError::InvalidOther(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMap(_) => write!(f, "invalid map"),
            Self::InvalidField(_) => write!(f, "invalid field"),
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
            Self::MissingId => write!(f, "missing ID"),
            Self::InvalidValues(_) => write!(f, "invalid values"),
            Self::InvalidOther(tag, _) => write!(f, "invalid other: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_other(src: &mut &[u8]) -> Result<(String, Map<Other>), ParseError> {
    const VALUES: &str = "Values";

    super::consume_prefix(src).map_err(ParseError::InvalidMap)?;

    let mut id = None;

    let mut other_fields = OtherFields::new();

    loop {
        let tag = parse_key(src)
            .map(Tag::from)
            .map_err(ParseError::InvalidKey)?;

        match tag {
            tag::ID => parse_id(src).and_then(|v| try_replace(&mut id, tag::ID, v))?,
            Tag::Other(t) => {
                if t.as_ref() == VALUES {
                    parse_values(src, &t)
                        .and_then(|value| try_insert(&mut other_fields, t, value))?;
                } else {
                    parse_other_value(src, &t)
                        .and_then(|value| try_insert(&mut other_fields, t, value))?;
                }
            }
        }

        if !super::field::consume_separator(src).map_err(ParseError::InvalidField)? {
            break;
        }
    }

    super::consume_suffix(src).map_err(ParseError::InvalidMap)?;

    let id = id.ok_or(ParseError::MissingId)?;

    Ok((
        id,
        Map {
            inner: Other,
            other_fields,
        },
    ))
}

fn parse_id(src: &mut &[u8]) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(ParseError::InvalidValue)
}

fn parse_values(
    src: &mut &[u8],
    tag: &map::tag::Other<tag::StandardTag>,
) -> Result<String, ParseError> {
    const PREFIX: u8 = b'[';
    const SUFFIX: u8 = b']';

    let is_delimited = src.first().map(|&b| b == PREFIX).unwrap_or_default();

    if is_delimited {
        if let Some(i) = src.iter().position(|&b| b == SUFFIX) {
            let (buf, rest) = src.split_at(i + 1);
            let s = str::from_utf8(buf).map_err(ParseError::InvalidValues)?;
            *src = rest;
            return Ok(s.into());
        }
    }

    parse_other_value(src, tag)
}

fn parse_other_value(
    src: &mut &[u8],
    tag: &map::tag::Other<tag::StandardTag>,
) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::InvalidOther(tag.clone(), e))
}

fn try_replace<T>(option: &mut Option<T>, tag: Tag, value: T) -> Result<(), ParseError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(ParseError::DuplicateTag(tag))
    }
}

fn try_insert(
    other_fields: &mut OtherFields<tag::StandardTag>,
    tag: map::tag::Other<tag::StandardTag>,
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
            Err(ParseError::DuplicateTag(Tag::Other(t)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_other() {
        let mut src = &br#"<ID=noodles>"#[..];

        let id = String::from("noodles");
        let map = Map::<Other>::new();
        let expected = (id, map);

        assert_eq!(parse_other(&mut src), Ok(expected));
    }
}
