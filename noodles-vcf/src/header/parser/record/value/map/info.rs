use std::{error, fmt, num};

use crate::header::{
    number,
    record::value::{
        map::{
            self,
            info::{tag, ty, Tag, Type},
            Info, OtherFields,
        },
        Map,
    },
    FileFormat, Number,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseErrorKind {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    MissingId,
    MissingNumber,
    InvalidNumber(number::ParseError),
    MissingType,
    InvalidType(ty::ParseError),
    MissingDescription,
    InvalidIdx(num::ParseIntError),
    DuplicateTag(Tag),
}

/// An error returned when a VCF header record info map value fails to parse.
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
            ParseErrorKind::InvalidNumber(e) => Some(e),
            ParseErrorKind::InvalidType(e) => Some(e),
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
            ParseErrorKind::MissingNumber => write!(f, "missing number"),
            ParseErrorKind::InvalidNumber(_) => write!(f, "invalid number"),
            ParseErrorKind::MissingType => write!(f, "missing type"),
            ParseErrorKind::InvalidType(_) => write!(f, "invalid type"),
            ParseErrorKind::MissingDescription => write!(f, "missing description"),
            ParseErrorKind::InvalidIdx(_) => write!(f, "invalid IDX"),
            ParseErrorKind::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_info(
    src: &mut &[u8],
    _file_format: FileFormat,
) -> Result<(String, Map<Info>), ParseError> {
    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut number = None;
    let mut ty = None;
    let mut description = None;
    let mut idx = None;

    let mut other_fields = OtherFields::new();

    while let Some((raw_key, raw_value)) = super::split_field(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidField(e)))?
    {
        match Tag::from(raw_key) {
            tag::ID => try_replace(&mut id, &None, tag::ID, raw_value.into())?,
            tag::NUMBER => parse_number(&raw_value, &id)
                .and_then(|v| try_replace(&mut number, &id, tag::NUMBER, v))?,
            tag::TYPE => {
                parse_type(&raw_value, &id)
                    .and_then(|v| try_replace(&mut ty, &id, tag::TYPE, v))?;
            }
            tag::DESCRIPTION => {
                try_replace(&mut description, &id, tag::DESCRIPTION, raw_value.into())?;
            }
            tag::IDX => {
                parse_idx(&raw_value, &id).and_then(|v| try_replace(&mut idx, &id, tag::IDX, v))?;
            }
            Tag::Other(t) => try_insert(&mut other_fields, &id, t, raw_value.into())?,
        }
    }

    super::consume_suffix(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidMap(e)))?;

    let id = id.ok_or_else(|| ParseError::new(None, ParseErrorKind::MissingId))?;
    let number =
        number.ok_or_else(|| ParseError::new(Some(id.clone()), ParseErrorKind::MissingNumber))?;
    let ty = ty.ok_or_else(|| ParseError::new(Some(id.clone()), ParseErrorKind::MissingType))?;
    let description = description
        .ok_or_else(|| ParseError::new(Some(id.clone()), ParseErrorKind::MissingDescription))?;

    Ok((
        id,
        Map {
            inner: Info {
                number,
                ty,
                description,
                idx,
            },
            other_fields,
        },
    ))
}

fn parse_number(s: &str, id: &Option<String>) -> Result<Number, ParseError> {
    s.parse()
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidNumber(e)))
}

fn parse_type(s: &str, id: &Option<String>) -> Result<Type, ParseError> {
    s.parse()
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidType(e)))
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
            let (t, _) = entry.swap_remove_entry();
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
    use crate::variant::record_buf::info::field::key;

    #[test]
    fn test_parse_info() {
        let mut src =
            &br#"<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#[..];
        let file_format = FileFormat::new(4, 4);

        let id = key::SAMPLES_WITH_DATA_COUNT;
        let map = Map::<Info>::from(id);
        let expected = (String::from(id), map);

        assert_eq!(parse_info(&mut src, file_format), Ok(expected));
    }
}
