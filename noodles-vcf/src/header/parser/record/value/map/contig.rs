use std::{error, fmt, num};

use crate::header::record::value::{
    map::{
        self,
        contig::{tag, Tag},
        Contig, OtherFields,
    },
    Map,
};

#[derive(Clone, Debug, Eq, PartialEq)]
enum ParseErrorKind {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    MissingId,
    InvalidLength(num::ParseIntError),
    InvalidIdx(num::ParseIntError),
    DuplicateTag(Tag),
}

/// An error returned when a VCF header record contig map value fails to parse.
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
            ParseErrorKind::InvalidLength(e) => Some(e),
            ParseErrorKind::InvalidIdx(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ParseErrorKind::InvalidMap(_) => write!(f, "invalid map"),
            ParseErrorKind::InvalidField(_) => write!(f, "invalid field"),
            ParseErrorKind::MissingId => write!(f, "missing ID"),
            ParseErrorKind::InvalidLength(_) => write!(f, "invalid length"),
            ParseErrorKind::InvalidIdx(_) => write!(f, "invalid IDX"),
            ParseErrorKind::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_contig(src: &mut &[u8]) -> Result<(String, Map<Contig>), ParseError> {
    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut length = None;
    let mut md5 = None;
    let mut url = None;
    let mut idx = None;

    let mut other_fields = OtherFields::new();

    while let Some((raw_key, raw_value)) = super::split_field(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidField(e)))?
    {
        match Tag::from(raw_key) {
            tag::ID => try_replace(&mut id, &None, tag::ID, raw_value.into())?,
            tag::LENGTH => parse_length(&raw_value, &id)
                .and_then(|v| try_replace(&mut length, &id, tag::LENGTH, v))?,
            tag::MD5 => try_replace(&mut md5, &id, tag::MD5, raw_value.into())?,
            tag::URL => try_replace(&mut url, &id, tag::URL, raw_value.into())?,
            tag::IDX => {
                parse_idx(&raw_value, &id).and_then(|v| try_replace(&mut idx, &id, tag::IDX, v))?;
            }
            Tag::Other(t) => try_insert(&mut other_fields, &id, t, raw_value.into())?,
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

fn parse_length(s: &str, id: &Option<String>) -> Result<usize, ParseError> {
    s.parse()
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidLength(e)))
}

fn parse_idx(s: &str, id: &Option<String>) -> Result<usize, ParseError> {
    s.parse()
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidIdx(e)))
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
    fn test_parse_contig() {
        let mut src = &br#"<ID=sq0>"#[..];

        let id = String::from("sq0");
        let map = Map::<Contig>::new();
        let expected = (id, map);

        assert_eq!(parse_contig(&mut src), Ok(expected));
    }
}
