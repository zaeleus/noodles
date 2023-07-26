use std::{error, fmt};

use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            read_group::{platform, tag, Platform, Tag},
            OtherFields, ReadGroup,
        },
        Map,
    },
};

/// An error returned when a SAM header read group record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),

    /// The predicted median insert size is invalid.
    InvalidPredictedMedianInsertSize(lexical_core::Error),
    /// The platform is invalid.
    InvalidPlatform(platform::ParseError),

    MissingField(Tag),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidPredictedMedianInsertSize(e) => Some(e),
            Self::InvalidPlatform(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(_) => write!(f, "invalid field"),
            Self::InvalidPredictedMedianInsertSize(_) => {
                write!(f, "invalid predicted median insert size")
            }
            Self::InvalidPlatform(_) => write!(f, "invalid platform"),
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_read_group(
    src: &mut &[u8],
    ctx: &Context,
) -> Result<(String, Map<ReadGroup>), ParseError> {
    let mut id = None;
    let mut barcode = None;
    let mut sequencing_center = None;
    let mut description = None;
    let mut produced_at = None;
    let mut flow_order = None;
    let mut key_sequence = None;
    let mut library = None;
    let mut program = None;
    let mut predicted_median_insert_size = None;
    let mut platform = None;
    let mut platform_model = None;
    let mut platform_unit = None;
    let mut sample = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidField)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::ID => parse_string(src).and_then(|v| try_replace(&mut id, ctx, tag::ID, v))?,
            tag::BARCODE => {
                parse_string(src).and_then(|v| try_replace(&mut barcode, ctx, tag::BARCODE, v))?
            }
            tag::SEQUENCING_CENTER => parse_string(src).and_then(|v| {
                try_replace(&mut sequencing_center, ctx, tag::SEQUENCING_CENTER, v)
            })?,
            tag::DESCRIPTION => parse_string(src)
                .and_then(|v| try_replace(&mut description, ctx, tag::DESCRIPTION, v))?,
            tag::PRODUCED_AT => parse_string(src)
                .and_then(|v| try_replace(&mut produced_at, ctx, tag::PRODUCED_AT, v))?,
            tag::FLOW_ORDER => parse_string(src)
                .and_then(|v| try_replace(&mut flow_order, ctx, tag::FLOW_ORDER, v))?,
            tag::KEY_SEQUENCE => parse_string(src)
                .and_then(|v| try_replace(&mut key_sequence, ctx, tag::KEY_SEQUENCE, v))?,
            tag::LIBRARY => {
                parse_string(src).and_then(|v| try_replace(&mut library, ctx, tag::LIBRARY, v))?
            }
            tag::PROGRAM => {
                parse_string(src).and_then(|v| try_replace(&mut program, ctx, tag::PROGRAM, v))?
            }
            tag::PREDICTED_MEDIAN_INSERT_SIZE => {
                parse_predicted_median_insert_size(src).and_then(|v| {
                    try_replace(
                        &mut predicted_median_insert_size,
                        ctx,
                        tag::PREDICTED_MEDIAN_INSERT_SIZE,
                        v,
                    )
                })?
            }
            tag::PLATFORM => parse_platform(src)
                .and_then(|v| try_replace(&mut platform, ctx, tag::PLATFORM, v))?,
            tag::PLATFORM_MODEL => parse_string(src)
                .and_then(|v| try_replace(&mut platform_model, ctx, tag::PLATFORM_MODEL, v))?,
            tag::PLATFORM_UNIT => parse_string(src)
                .and_then(|v| try_replace(&mut platform_unit, ctx, tag::PLATFORM_UNIT, v))?,
            tag::SAMPLE => {
                parse_string(src).and_then(|v| try_replace(&mut sample, ctx, tag::SAMPLE, v))?
            }
            Tag::Other(t) => parse_value(src)
                .map_err(ParseError::InvalidField)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let id = id.ok_or(ParseError::MissingField(tag::ID))?;

    Ok((
        id,
        Map {
            inner: ReadGroup {
                barcode,
                sequencing_center,
                description,
                produced_at,
                flow_order,
                key_sequence,
                library,
                program,
                predicted_median_insert_size,
                platform,
                platform_model,
                platform_unit,
                sample,
            },
            other_fields,
        },
    ))
}

fn parse_string(src: &mut &[u8]) -> Result<String, ParseError> {
    parse_value(src)
        .map(String::from)
        .map_err(ParseError::InvalidField)
}

fn parse_predicted_median_insert_size(src: &mut &[u8]) -> Result<i32, ParseError> {
    let (n, i) =
        lexical_core::parse_partial(src).map_err(ParseError::InvalidPredictedMedianInsertSize)?;

    *src = &src[i..];

    Ok(n)
}

fn parse_platform(src: &mut &[u8]) -> Result<Platform, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidField)
        .and_then(|s| s.parse().map_err(ParseError::InvalidPlatform))
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
    fn test_parse_read_group() {
        let mut src = &b"\tID:rg0"[..];
        let ctx = Context::default();
        let actual = parse_read_group(&mut src, &ctx);
        let expected = (String::from("rg0"), Map::<ReadGroup>::default());
        assert_eq!(actual, Ok(expected));
    }

    #[test]
    fn test_parse_read_group_with_missing_id() {
        let mut src = &b"\tPL:ILLUMINA"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_read_group(&mut src, &ctx),
            Err(ParseError::MissingField(tag::ID))
        );
    }
}
