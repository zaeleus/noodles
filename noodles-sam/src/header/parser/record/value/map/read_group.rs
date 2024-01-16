use std::{error, fmt};

use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value, value};
use crate::header::{
    parser::Context,
    record::value::{
        map::{
            self,
            read_group::{platform, tag, Platform, Tag},
            tag::Other,
            OtherFields, ReadGroup,
        },
        Map,
    },
};

/// An error returned when a SAM header read group record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),
    InvalidTag(super::field::tag::ParseError),
    InvalidValue(value::ParseError),
    MissingId,
    InvalidId(value::ParseError),
    InvalidBarcode(value::ParseError),
    InvalidSequencingCenter(value::ParseError),
    InvalidDescription(value::ParseError),
    InvalidProducedAt(value::ParseError),
    InvalidFlowOrder(value::ParseError),
    InvalidKeySequence(value::ParseError),
    InvalidLibrary(value::ParseError),
    InvalidProgram(value::ParseError),
    InvalidPredictedMedianInsertSize(lexical_core::Error),
    InvalidPlatform(platform::ParseError),
    InvalidPlatformModel(value::ParseError),
    InvalidPlatformUnit(value::ParseError),
    InvalidSample(value::ParseError),
    InvalidOther(Other<tag::Standard>, value::ParseError),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidTag(e) => Some(e),
            Self::InvalidId(e) => Some(e),
            Self::InvalidBarcode(e) => Some(e),
            Self::InvalidSequencingCenter(e) => Some(e),
            Self::InvalidDescription(e) => Some(e),
            Self::InvalidProducedAt(e) => Some(e),
            Self::InvalidFlowOrder(e) => Some(e),
            Self::InvalidKeySequence(e) => Some(e),
            Self::InvalidLibrary(e) => Some(e),
            Self::InvalidProgram(e) => Some(e),
            Self::InvalidPredictedMedianInsertSize(e) => Some(e),
            Self::InvalidPlatform(e) => Some(e),
            Self::InvalidPlatformModel(e) => Some(e),
            Self::InvalidPlatformUnit(e) => Some(e),
            Self::InvalidSample(e) => Some(e),
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
            Self::MissingId => write!(f, "missing ID"),
            Self::InvalidId(_) => write!(f, "invalid ID"),
            Self::InvalidBarcode(_) => write!(f, "invalid barcode ({})", tag::BARCODE),
            Self::InvalidSequencingCenter(_) => {
                write!(f, "invalid sequencing center ({})", tag::SEQUENCING_CENTER)
            }
            Self::InvalidDescription(_) => write!(f, "invalid description ({})", tag::DESCRIPTION),
            Self::InvalidProducedAt(_) => write!(f, "invalid produced at ({})", tag::PRODUCED_AT),
            Self::InvalidFlowOrder(_) => write!(f, "invalid flow order ({})", tag::FLOW_ORDER),
            Self::InvalidKeySequence(_) => {
                write!(f, "invalid key sequence ({})", tag::KEY_SEQUENCE)
            }
            Self::InvalidLibrary(_) => write!(f, "invalid library ({})", tag::LIBRARY),
            Self::InvalidProgram(_) => write!(f, "invalid program ({})", tag::PROGRAM),
            Self::InvalidPredictedMedianInsertSize(_) => {
                write!(
                    f,
                    "invalid predicted median insert size ({})",
                    tag::PREDICTED_MEDIAN_INSERT_SIZE
                )
            }
            Self::InvalidPlatform(_) => write!(f, "invalid platform ({})", tag::PLATFORM),
            Self::InvalidPlatformModel(_) => {
                write!(f, "invalid platform model ({})", tag::PLATFORM_MODEL)
            }
            Self::InvalidPlatformUnit(_) => {
                write!(f, "invalid platform unit ({})", tag::PLATFORM_UNIT)
            }
            Self::InvalidSample(_) => write!(f, "invalid sample ({})", tag::SAMPLE),
            Self::InvalidOther(tag, _) => write!(f, "invalid other ({tag})"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_read_group(
    src: &mut &[u8],
    ctx: &Context,
) -> Result<(Vec<u8>, Map<ReadGroup>), ParseError> {
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
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::ID => parse_id(src).and_then(|v| try_replace(&mut id, ctx, tag::ID, v))?,
            tag::BARCODE => {
                parse_barcode(src).and_then(|v| try_replace(&mut barcode, ctx, tag::BARCODE, v))?;
            }
            tag::SEQUENCING_CENTER => parse_sequencing_center(src).and_then(|v| {
                try_replace(&mut sequencing_center, ctx, tag::SEQUENCING_CENTER, v)
            })?,
            tag::DESCRIPTION => parse_description(src)
                .and_then(|v| try_replace(&mut description, ctx, tag::DESCRIPTION, v))?,
            tag::PRODUCED_AT => parse_produced_at(src)
                .and_then(|v| try_replace(&mut produced_at, ctx, tag::PRODUCED_AT, v))?,
            tag::FLOW_ORDER => parse_flow_order(src)
                .and_then(|v| try_replace(&mut flow_order, ctx, tag::FLOW_ORDER, v))?,
            tag::KEY_SEQUENCE => parse_key_sequence(src)
                .and_then(|v| try_replace(&mut key_sequence, ctx, tag::KEY_SEQUENCE, v))?,
            tag::LIBRARY => {
                parse_library(src).and_then(|v| try_replace(&mut library, ctx, tag::LIBRARY, v))?;
            }
            tag::PROGRAM => {
                parse_program(src).and_then(|v| try_replace(&mut program, ctx, tag::PROGRAM, v))?;
            }
            tag::PREDICTED_MEDIAN_INSERT_SIZE => {
                parse_predicted_median_insert_size(src).and_then(|v| {
                    try_replace(
                        &mut predicted_median_insert_size,
                        ctx,
                        tag::PREDICTED_MEDIAN_INSERT_SIZE,
                        v,
                    )
                })?;
            }
            tag::PLATFORM => parse_platform(src)
                .and_then(|v| try_replace(&mut platform, ctx, tag::PLATFORM, v))?,
            tag::PLATFORM_MODEL => parse_platform_model(src)
                .and_then(|v| try_replace(&mut platform_model, ctx, tag::PLATFORM_MODEL, v))?,
            tag::PLATFORM_UNIT => parse_platform_unit(src)
                .and_then(|v| try_replace(&mut platform_unit, ctx, tag::PLATFORM_UNIT, v))?,
            tag::SAMPLE => {
                parse_sample(src).and_then(|v| try_replace(&mut sample, ctx, tag::SAMPLE, v))?;
            }
            Tag::Other(t) => parse_other(src, t)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let id = id.ok_or(ParseError::MissingId)?;

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

fn parse_id(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidId)
}

fn parse_barcode(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidBarcode)
}

fn parse_sequencing_center(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidSequencingCenter)
}

fn parse_description(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidDescription)
}

fn parse_produced_at(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidProducedAt)
}

fn parse_flow_order(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidFlowOrder)
}

fn parse_key_sequence(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidKeySequence)
}

fn parse_library(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidLibrary)
}

fn parse_program(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidProgram)
}

fn parse_predicted_median_insert_size(src: &mut &[u8]) -> Result<i32, ParseError> {
    let (n, i) =
        lexical_core::parse_partial(src).map_err(ParseError::InvalidPredictedMedianInsertSize)?;

    *src = &src[i..];

    Ok(n)
}

fn parse_platform(src: &mut &[u8]) -> Result<Platform, ParseError> {
    parse_value(src)
        .map_err(ParseError::InvalidValue)
        .and_then(|s| s.parse().map_err(ParseError::InvalidPlatform))
}

fn parse_platform_model(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidPlatformModel)
}

fn parse_platform_unit(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidPlatformUnit)
}

fn parse_sample(src: &mut &[u8]) -> Result<Vec<u8>, ParseError> {
    parse_value(src)
        .map(Vec::from)
        .map_err(ParseError::InvalidSample)
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
    use super::*;

    #[test]
    fn test_parse_read_group() {
        let mut src = &b"\tID:rg0"[..];
        let ctx = Context::default();
        let actual = parse_read_group(&mut src, &ctx);
        let expected = (Vec::from("rg0"), Map::<ReadGroup>::default());
        assert_eq!(actual, Ok(expected));
    }

    #[test]
    fn test_parse_read_group_with_missing_id() {
        let mut src = &b"\tPL:ILLUMINA"[..];
        let ctx = Context::default();
        assert_eq!(parse_read_group(&mut src, &ctx), Err(ParseError::MissingId));
    }
}
