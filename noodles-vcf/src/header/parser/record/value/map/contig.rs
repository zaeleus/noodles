use std::{error, fmt, num};

use super::field::{parse_key, parse_value};
use crate::header::record::value::{
    map::{
        self,
        contig::{
            name::{self, Name},
            tag, Tag,
        },
        Contig, OtherFields,
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
    InvalidId(name::ParseError),
    InvalidLength(num::ParseIntError),
    InvalidIdx(num::ParseIntError),
    InvalidOther(
        map::tag::Other<tag::Standard>,
        super::field::value::ParseError,
    ),
    DuplicateTag(Tag),
}

/// An error returned when a VCF header record contig map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError {
    id: Option<Name>,
    kind: ParseErrorKind,
}

impl ParseError {
    fn new(id: Option<Name>, kind: ParseErrorKind) -> Self {
        Self { id, kind }
    }

    pub(crate) fn id(&self) -> Option<&Name> {
        self.id.as_ref()
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match &self.kind {
            ParseErrorKind::InvalidMap(e) => Some(e),
            ParseErrorKind::InvalidField(e) => Some(e),
            ParseErrorKind::InvalidKey(e) => Some(e),
            ParseErrorKind::InvalidValue(e) => Some(e),
            ParseErrorKind::InvalidId(e) => Some(e),
            ParseErrorKind::InvalidLength(e) => Some(e),
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
            ParseErrorKind::InvalidId(_) => write!(f, "invalid ID"),
            ParseErrorKind::InvalidLength(_) => write!(f, "invalid length"),
            ParseErrorKind::InvalidIdx(_) => write!(f, "invalid IDX"),
            ParseErrorKind::InvalidOther(tag, _) => write!(f, "invalid other: {tag}"),
            ParseErrorKind::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_contig(src: &mut &[u8]) -> Result<(Name, Map<Contig>), ParseError> {
    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut length = None;
    let mut md5 = None;
    let mut url = None;
    let mut idx = None;

    let mut other_fields = OtherFields::new();

    loop {
        let tag = parse_key(src)
            .map(Tag::from)
            .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidKey(e)))?;

        match tag {
            tag::ID => parse_id(src, &id).and_then(|v| try_replace(&mut id, &None, tag::ID, v))?,
            tag::LENGTH => parse_length(src, &id)
                .and_then(|v| try_replace(&mut length, &id, tag::LENGTH, v))?,
            tag::MD5 => {
                parse_md5(src, &id).and_then(|v| try_replace(&mut md5, &id, tag::MD5, v))?
            }
            tag::URL => {
                parse_url(src, &id).and_then(|v| try_replace(&mut url, &id, tag::URL, v))?
            }
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

    Ok((
        id,
        Map {
            inner: Contig {
                length,
                md5,
                url,
                idx,
            },
            other_fields,
        },
    ))
}

fn parse_id(src: &mut &[u8], id: &Option<Name>) -> Result<Name, ParseError> {
    parse_value(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidId(e)))
        })
}

fn parse_length(src: &mut &[u8], id: &Option<Name>) -> Result<usize, ParseError> {
    parse_value(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidLength(e)))
        })
}

fn parse_md5(src: &mut &[u8], id: &Option<Name>) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
}

fn parse_url(src: &mut &[u8], id: &Option<Name>) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
}

fn parse_idx(src: &mut &[u8], id: &Option<Name>) -> Result<usize, ParseError> {
    parse_value(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidIdx(e)))
        })
}

fn parse_other(
    src: &mut &[u8],
    id: &Option<Name>,
    tag: &map::tag::Other<tag::Standard>,
) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidOther(tag.clone(), e)))
}

fn try_replace<T>(
    option: &mut Option<T>,
    id: &Option<Name>,
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
    id: &Option<Name>,
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
    fn test_parse_contig() -> Result<(), name::ParseError> {
        let mut src = &br#"<ID=sq0>"#[..];

        let id = "sq0".parse()?;
        let map = Map::<Contig>::new();
        let expected = (id, map);

        assert_eq!(parse_contig(&mut src), Ok(expected));

        Ok(())
    }
}
