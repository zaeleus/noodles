use std::{error, fmt, str};

use super::field::{parse_key, parse_value};
use crate::header::{
    record::value::{
        map::{
            self,
            other::{tag, Tag},
            Other, OtherFields,
        },
        Map,
    },
    FileFormat,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseErrorKind {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    InvalidKey(super::field::key::ParseError),
    InvalidValue(super::field::value::ParseError),
    MissingId,
    InvalidValues(str::Utf8Error),
    InvalidOther(
        map::tag::Other<tag::Standard>,
        super::field::value::ParseError,
    ),
    DuplicateTag(Tag),
}

/// An error returned when a VCF header record other map value fails to parse.
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
            ParseErrorKind::InvalidValues(e) => Some(e),
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
            ParseErrorKind::InvalidValues(_) => write!(f, "invalid values"),
            ParseErrorKind::InvalidOther(tag, _) => write!(f, "invalid other: {tag}"),
            ParseErrorKind::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_other(src: &mut &[u8]) -> Result<(String, Map<Other>), ParseError> {
    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut other_fields = OtherFields::new();

    while let Some((raw_key, raw_value)) = super::split_field(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidField(e)))?
    {
        match Tag::from(raw_key) {
            tag::ID => try_replace(&mut id, &None, tag::ID, raw_value.into())?,
            Tag::Other(t) => try_insert(&mut other_fields, &id, t, raw_value.into())?,
        }
    }

    super::consume_suffix(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidMap(e)))?;

    let id = id.ok_or_else(|| ParseError::new(None, ParseErrorKind::MissingId))?;

    Ok((
        id,
        Map {
            inner: Other::default(),
            other_fields,
        },
    ))
}

pub fn parse_meta(
    src: &mut &[u8],
    file_format: FileFormat,
) -> Result<(String, Map<Other>), ParseError> {
    const VALUES: &str = "Values";
    const VCF_4_3: FileFormat = FileFormat::new(4, 3);

    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut other_fields = OtherFields::new();

    loop {
        let tag = parse_key(src)
            .map(Tag::from)
            .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidKey(e)))?;

        match tag {
            tag::ID => parse_id(src, &id).and_then(|v| try_replace(&mut id, &None, tag::ID, v))?,
            Tag::Other(t) => {
                if file_format >= VCF_4_3 && t == VALUES {
                    parse_values(src, &id, &t)
                        .and_then(|value| try_insert(&mut other_fields, &id, t, value))?;
                } else {
                    parse_other_value(src, &id, &t)
                        .and_then(|value| try_insert(&mut other_fields, &id, t, value))?;
                }
            }
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
            inner: Other::default(),
            other_fields,
        },
    ))
}

pub fn parse_pedigree(
    src: &mut &[u8],
    file_format: FileFormat,
) -> Result<(String, Map<Other>), ParseError> {
    const CHILD: &str = "Child";
    const DERIVED: &str = "Derived";
    const VCF_4_3: FileFormat = FileFormat::new(4, 3);

    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id_tag = tag::ID;
    let mut id = None;

    let mut other_fields = OtherFields::new();

    loop {
        let tag = parse_key(src)
            .map(Tag::from)
            .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidKey(e)))?;

        match tag {
            tag::ID => parse_id(src, &id).and_then(|v| try_replace(&mut id, &None, tag::ID, v))?,
            Tag::Other(t) => {
                if file_format < VCF_4_3 && matches!(t.as_ref(), CHILD | DERIVED) {
                    id_tag = Tag::Other(t);
                    parse_id(src, &id).and_then(|v| try_replace(&mut id, &None, tag::ID, v))?;
                } else {
                    parse_other_value(src, &id, &t)
                        .and_then(|value| try_insert(&mut other_fields, &id, t, value))?;
                }
            }
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
            inner: Other { id_tag },
            other_fields,
        },
    ))
}

fn parse_id(src: &mut &[u8], id: &Option<String>) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValue(e)))
}

fn parse_values(
    src: &mut &[u8],
    id: &Option<String>,
    tag: &map::tag::Other<tag::Standard>,
) -> Result<String, ParseError> {
    use memchr::memchr;

    const PREFIX: u8 = b'[';
    const SUFFIX: u8 = b']';

    let is_delimited = src.first().map(|&b| b == PREFIX).unwrap_or_default();

    if is_delimited {
        if let Some(i) = memchr(SUFFIX, src) {
            let (buf, rest) = src.split_at(i + 1);

            let s = str::from_utf8(buf)
                .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidValues(e)))?;

            *src = rest;

            return Ok(s.into());
        }
    }

    parse_other_value(src, id, tag)
}

fn parse_other_value(
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
        Entry::Occupied(entry) => Err(ParseError::new(
            id.clone(),
            ParseErrorKind::DuplicateTag(Tag::Other(entry.key().clone())),
        )),
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

    #[test]
    fn test_parse_meta() -> Result<(), Box<dyn std::error::Error>> {
        const VCF_4_2: FileFormat = FileFormat::new(4, 2);
        const VCF_4_3: FileFormat = FileFormat::new(4, 3);

        let mut src = &b"<ID=Assay,Values=[WholeGenome, Exome]>"[..];
        assert_eq!(
            parse_meta(&mut src, VCF_4_3),
            Ok((
                String::from("Assay"),
                Map::<Other>::builder()
                    .insert("Values".parse()?, "[WholeGenome, Exome]")
                    .build()?
            ))
        );

        let mut src = &b"<ID=Assay,Values=[WholeGenome, Exome]>"[..];
        assert!(matches!(
            parse_meta(&mut src, VCF_4_2),
            Err(ParseError {
                id,
                kind: ParseErrorKind::InvalidKey(_)
            }) if id == Some(String::from("Assay"))
        ));

        Ok(())
    }

    #[test]
    fn test_parse_pedigree() -> Result<(), Box<dyn std::error::Error>> {
        const VCF_4_2: FileFormat = FileFormat::new(4, 2);
        const VCF_4_3: FileFormat = FileFormat::new(4, 3);

        let mut src = &b"<ID=CID>"[..];
        assert_eq!(
            parse_pedigree(&mut src, VCF_4_3),
            Ok((String::from("CID"), Map::<Other>::new()))
        );

        let mut src = &b"<Derived=DID,Original=OID>"[..];
        assert_eq!(
            parse_pedigree(&mut src, VCF_4_3),
            Err(ParseError::new(None, ParseErrorKind::MissingId))
        );

        let mut src = &b"<Child=CID,Mother=MID,Father=FID>"[..];
        assert_eq!(
            parse_pedigree(&mut src, VCF_4_3),
            Err(ParseError::new(None, ParseErrorKind::MissingId))
        );

        let mut src = &b"<Derived=DID,Original=OID>"[..];
        let (actual_id, actual_map) = parse_pedigree(&mut src, VCF_4_2)?;

        let expected_id = String::from("DID");
        let mut expected_map = Map::<Other>::builder()
            .insert("Original".parse()?, "OID")
            .build()?;
        expected_map.inner.id_tag = Tag::from("Derived");

        assert_eq!((actual_id, actual_map), (expected_id, expected_map));

        let mut src = &b"<Child=CID,Mother=MID,Father=FID>"[..];
        let (actual_id, actual_map) = parse_pedigree(&mut src, VCF_4_2)?;

        let expected_id = String::from("CID");
        let mut expected_map = Map::<Other>::builder()
            .insert("Mother".parse()?, "MID")
            .insert("Father".parse()?, "FID")
            .build()?;
        expected_map.inner.id_tag = Tag::from("Child");

        assert_eq!((actual_id, actual_map), (expected_id, expected_map));

        Ok(())
    }
}
