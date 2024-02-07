//! VCF record reader.

mod alternate_bases;
mod chromosome;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_bases;
mod samples;
mod value;

use std::{error, fmt};

use self::{
    alternate_bases::parse_alternate_bases, chromosome::parse_chromosome, filters::parse_filters,
    ids::parse_ids, info::parse_info, position::parse_position, quality_score::parse_quality_score,
    reference_bases::parse_reference_bases, samples::parse_samples,
};
use crate::{variant::RecordBuf, Header};

const MISSING: &str = ".";

/// An error when a raw VCF record fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The position is invalid.
    InvalidPosition(position::ParseError),
    /// The IDs are invalid.
    InvalidIds(ids::ParseError),
    /// The reference bases are invalid.
    InvalidReferenceBases(reference_bases::ParseError),
    /// The alternate bases are invalid.
    InvalidAlternateBases(alternate_bases::ParseError),
    /// The quality score is invalid.
    InvalidQualityScore(quality_score::ParseError),
    /// The filters are invalid.
    InvalidFilters(filters::ParseError),
    /// The info is invalid.
    InvalidInfo(info::ParseError),
    /// The samples are invalid.
    InvalidSamples(samples::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidPosition(e) => Some(e),
            Self::InvalidIds(e) => Some(e),
            Self::InvalidReferenceBases(e) => Some(e),
            Self::InvalidAlternateBases(e) => Some(e),
            Self::InvalidQualityScore(e) => Some(e),
            Self::InvalidFilters(e) => Some(e),
            Self::InvalidInfo(e) => Some(e),
            Self::InvalidSamples(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPosition(_) => write!(f, "invalid position"),
            Self::InvalidIds(_) => write!(f, "invalid IDs"),
            Self::InvalidReferenceBases(_) => write!(f, "invalid reference bases"),
            Self::InvalidAlternateBases(_) => write!(f, "invalid alternate bases"),
            Self::InvalidQualityScore(_) => write!(f, "invalid quality score"),
            Self::InvalidFilters(_) => write!(f, "invalid filters"),
            Self::InvalidInfo(_) => write!(f, "invalid info"),
            Self::InvalidSamples(_) => write!(f, "invalid samples"),
        }
    }
}

pub(crate) fn parse_record(
    mut s: &str,
    header: &Header,
    record: &mut RecordBuf,
) -> Result<(), ParseError> {
    let field = next_field(&mut s);
    parse_chromosome(field, record.chromosome_mut());

    let field = next_field(&mut s);
    *record.position_mut() = parse_position(field).map_err(ParseError::InvalidPosition)?;

    record.ids_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_ids(field, record.ids_mut()).map_err(ParseError::InvalidIds)?;
    }

    let field = next_field(&mut s);
    parse_reference_bases(field, record.reference_bases_mut())
        .map_err(ParseError::InvalidReferenceBases)?;

    record.alternate_bases_mut().as_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_alternate_bases(field, record.alternate_bases_mut())
            .map_err(ParseError::InvalidAlternateBases)?;
    }

    let field = next_field(&mut s);
    *record.quality_score_mut() = match field {
        MISSING => None,
        _ => parse_quality_score(field)
            .map(Some)
            .map_err(ParseError::InvalidQualityScore)?,
    };

    let field = next_field(&mut s);
    match field {
        MISSING => {
            record.filters_mut().take();
        }
        _ => parse_filters(field, record.filters_mut()).map_err(ParseError::InvalidFilters)?,
    }

    record.info_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_info(header, field, record.info_mut()).map_err(ParseError::InvalidInfo)?;
    }

    parse_samples(header, s, record.samples_mut()).map_err(ParseError::InvalidSamples)?;

    Ok(())
}

fn next_field<'a>(s: &mut &'a str) -> &'a str {
    const DELIMITER: char = '\t';

    let (field, rest) = s
        .split_once(DELIMITER)
        .unwrap_or_else(|| s.split_at(s.len()));

    *s = rest;

    field
}
