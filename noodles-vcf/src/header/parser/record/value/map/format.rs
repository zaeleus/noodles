use std::{error, fmt, num};

use super::field::{parse_key, parse_value};
use crate::{
    header::{
        number,
        record::value::{
            map::{
                self,
                format::{tag, ty, Tag, Type},
                Format, OtherFields,
            },
            Map,
        },
        Number,
    },
    record::genotypes::keys::{key, Key},
};

/// An error returned when a VCF header record format map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    InvalidKey(super::field::key::ParseError),
    InvalidValue(super::field::value::ParseError),
    MissingId,
    InvalidId(key::ParseError),
    MissingNumber,
    InvalidNumber(number::ParseError),
    MissingType,
    InvalidType(ty::ParseError),
    MissingDescription,
    InvalidDescription,
    InvalidIdx(num::ParseIntError),
    InvalidOther(
        map::tag::Other<tag::StandardTag>,
        super::field::value::ParseError,
    ),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidMap(e) => Some(e),
            Self::InvalidField(e) => Some(e),
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValue(e) => Some(e),
            Self::InvalidId(e) => Some(e),
            Self::InvalidNumber(e) => Some(e),
            Self::InvalidType(e) => Some(e),
            Self::InvalidIdx(e) => Some(e),
            Self::InvalidOther(_, e) => Some(e),
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
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::MissingNumber => write!(f, "missing number"),
            Self::InvalidNumber(_) => write!(f, "invalid number"),
            Self::MissingType => write!(f, "missing type"),
            Self::InvalidType(_) => write!(f, "invalid type"),
            Self::MissingDescription => write!(f, "missing description"),
            Self::InvalidDescription => write!(f, "invalid description"),
            Self::InvalidIdx(_) => write!(f, "invalid IDX"),
            Self::InvalidOther(tag, _) => write!(f, "invalid other: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_format(src: &mut &[u8]) -> Result<(Key, Map<Format>), ParseError> {
    super::consume_prefix(src).map_err(ParseError::InvalidMap)?;

    let mut id = None;
    let mut number = None;
    let mut ty = None;
    let mut description = None;
    let mut idx = None;

    let mut other_fields = OtherFields::new();

    loop {
        let tag = parse_key(src)
            .map(Tag::from)
            .map_err(ParseError::InvalidKey)?;

        match tag {
            tag::ID => parse_id(src).and_then(|v| try_replace(&mut id, tag::ID, v))?,
            tag::NUMBER => {
                parse_number(src).and_then(|v| try_replace(&mut number, tag::NUMBER, v))?
            }
            tag::TYPE => parse_type(src).and_then(|v| try_replace(&mut ty, tag::TYPE, v))?,
            tag::DESCRIPTION => parse_description(src)
                .and_then(|v| try_replace(&mut description, tag::DESCRIPTION, v))?,
            tag::IDX => parse_idx(src).and_then(|v| try_replace(&mut idx, tag::IDX, v))?,
            Tag::Other(t) => {
                parse_other(src, &t).and_then(|value| try_insert(&mut other_fields, t, value))?
            }
        }

        if !super::field::consume_separator(src).map_err(ParseError::InvalidField)? {
            break;
        }
    }

    super::consume_suffix(src).map_err(ParseError::InvalidMap)?;

    let id = id.ok_or(ParseError::MissingId)?;
    let number = number.ok_or(ParseError::MissingNumber)?;
    let ty = ty.ok_or(ParseError::MissingType)?;
    let description = description.ok_or(ParseError::MissingDescription)?;

    Ok((
        id,
        Map {
            inner: Format {
                number,
                ty,
                description,
                idx,
            },
            other_fields,
        },
    ))
}

fn parse_id(src: &mut &[u8]) -> Result<Key, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidId))
}

fn parse_number(src: &mut &[u8]) -> Result<Number, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidNumber))
}

fn parse_type(src: &mut &[u8]) -> Result<Type, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidType))
}

fn parse_description(src: &mut &[u8]) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(ParseError::InvalidValue)
}

fn parse_idx(src: &mut &[u8]) -> Result<usize, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidIdx))
}

fn parse_other(
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
    fn test_parse_format() {
        let mut src = &br#"<ID=GT,Number=1,Type=String,Description="Genotype">"#[..];

        let id = key::GENOTYPE;
        let map = Map::<Format>::from(&id);
        let expected = (id, map);

        assert_eq!(parse_format(&mut src), Ok(expected));
    }
}
