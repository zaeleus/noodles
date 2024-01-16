mod length;

use std::{error, fmt};

use self::length::parse_length;
use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value, value};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            reference_sequence::{tag, Tag},
            tag::Other,
            OtherFields, ReferenceSequence,
        },
        Map,
    },
};

/// An error returned when a SAM header reference sequence record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),
    InvalidTag(super::field::tag::ParseError),
    InvalidValue(value::ParseError),
    MissingName,
    MissingLength,
    InvalidLength(length::ParseError),
    InvalidOther(Other<tag::Standard>, value::ParseError),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidTag(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
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
            Self::MissingName => write!(f, "missing name ({})", tag::NAME),
            Self::MissingLength => write!(f, "missing length ({})", tag::LENGTH),
            Self::InvalidLength(_) => write!(f, "invalid length ({})", tag::LENGTH),
            Self::InvalidOther(tag, _) => write!(f, "invalid other ({tag})"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_reference_sequence(
    src: &mut &[u8],
    ctx: &Context,
) -> Result<(Vec<u8>, Map<ReferenceSequence>), ParseError> {
    let mut name = None;
    let mut length = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::NAME => parse_name(src).and_then(|v| try_replace(&mut name, ctx, tag::NAME, v))?,
            tag::LENGTH => parse_length(src)
                .map_err(ParseError::InvalidLength)
                .and_then(|v| try_replace(&mut length, ctx, tag::LENGTH, v))?,
            Tag::Other(t) => parse_other(src, t)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let name = name.ok_or(ParseError::MissingName)?;
    let length = length.ok_or(ParseError::MissingLength)?;

    Ok((
        name,
        Map {
            inner: ReferenceSequence { length },
            other_fields,
        },
    ))
}

fn parse_name(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidValue)
}

fn parse_other(src: &mut &[u8], tag: Other<tag::Standard>) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
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
    V: Into<Vec<u8>>,
{
    if other_fields.insert(tag, value.into()).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(Tag::Other(tag)))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;

    #[test]
    fn test_parse_header() -> Result<(), map::reference_sequence::name::ParseError> {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let mut src = &b"\tSN:sq0\tLN:8"[..];
        let ctx = Context::default();
        let actual = parse_reference_sequence(&mut src, &ctx);

        let expected = (Vec::from("sq0"), Map::<ReferenceSequence>::new(SQ0_LN));

        assert_eq!(actual, Ok(expected));

        Ok(())
    }

    #[test]
    fn test_parse_header_with_missing_name() {
        let mut src = &b"\tLN:8"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_reference_sequence(&mut src, &ctx),
            Err(ParseError::MissingName)
        );
    }

    #[test]
    fn test_parse_header_with_missing_length() {
        let mut src = &b"\tSN:sq0"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_reference_sequence(&mut src, &ctx),
            Err(ParseError::MissingLength)
        );
    }

    #[test]
    fn test_parse_header_with_invalid_length() {
        let ctx = Context::default();

        let mut src = &b"\tSN:sq0\tLN:NA"[..];
        assert!(matches!(
            parse_reference_sequence(&mut src, &ctx),
            Err(ParseError::InvalidLength(_))
        ));

        let mut src = &b"\tSN:sq0\tLN:0"[..];
        assert!(matches!(
            parse_reference_sequence(&mut src, &ctx),
            Err(ParseError::InvalidLength(_))
        ));
    }
}
