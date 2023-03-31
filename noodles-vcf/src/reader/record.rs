mod chromosome;
mod filters;
mod genotypes;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_bases;

use std::{error, fmt, io};

use noodles_core as core;

use self::{
    chromosome::parse_chromosome, filters::parse_filters, genotypes::parse_genotypes,
    ids::parse_ids, info::parse_info, position::parse_position, quality_score::parse_quality_score,
    reference_bases::parse_reference_bases,
};
use crate::{record::AlternateBases, Header, Record};

const MISSING: &str = ".";

/// An error when a raw VCF record fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The chromosome is invalid.
    InvalidChromosome(chromosome::ParseError),
    /// The position is invalid.
    InvalidPosition(position::ParseError),
    /// The IDs are invalid.
    InvalidIds(ids::ParseError),
    /// The reference bases are invalid.
    InvalidReferenceBases(reference_bases::ParseError),
    /// The quality score is invalid.
    InvalidQualityScore(quality_score::ParseError),
    /// The filters are invalid.
    InvalidFilters(filters::ParseError),
    /// The info is invalid.
    InvalidInfo(info::ParseError),
    /// The genotypes are invalid.
    InvalidGenotypes(genotypes::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidChromosome(e) => Some(e),
            Self::InvalidPosition(e) => Some(e),
            Self::InvalidIds(e) => Some(e),
            Self::InvalidReferenceBases(e) => Some(e),
            Self::InvalidQualityScore(e) => Some(e),
            Self::InvalidFilters(e) => Some(e),
            Self::InvalidInfo(e) => Some(e),
            Self::InvalidGenotypes(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidChromosome(_) => write!(f, "invalid chromosome"),
            Self::InvalidPosition(_) => write!(f, "invalid position"),
            Self::InvalidIds(_) => write!(f, "invalid IDs"),
            Self::InvalidReferenceBases(_) => write!(f, "invalid reference bases"),
            Self::InvalidQualityScore(_) => write!(f, "invalid quality score"),
            Self::InvalidFilters(_) => write!(f, "invalid filters"),
            Self::InvalidInfo(_) => write!(f, "invalid info"),
            Self::InvalidGenotypes(_) => write!(f, "invalid genotypes"),
        }
    }
}

impl From<ParseError> for core::Error {
    fn from(e: ParseError) -> Self {
        Self::new(core::error::Kind::Parse, e)
    }
}

pub(super) fn parse_record(mut s: &str, header: &Header, record: &mut Record) -> io::Result<()> {
    let field = next_field(&mut s);
    parse_chromosome(field, record.chromosome_mut())
        .map_err(ParseError::InvalidChromosome)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let field = next_field(&mut s);
    *record.position_mut() = parse_position(field)
        .map_err(ParseError::InvalidPosition)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    record.ids_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_ids(field, record.ids_mut())
            .map_err(ParseError::InvalidIds)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    let field = next_field(&mut s);
    parse_reference_bases(field, record.reference_bases_mut())
        .map_err(ParseError::InvalidReferenceBases)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let field = next_field(&mut s);
    *record.alternate_bases_mut() = parse_alternate_bases(field)?;

    let field = next_field(&mut s);
    *record.quality_score_mut() = match field {
        MISSING => None,
        _ => parse_quality_score(field)
            .map(Some)
            .map_err(ParseError::InvalidQualityScore)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
    };

    let field = next_field(&mut s);
    match field {
        MISSING => {
            record.filters_mut().take();
        }
        _ => parse_filters(field, record.filters_mut())
            .map_err(ParseError::InvalidFilters)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
    }

    record.info_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_info(header, field, record.info_mut())
            .map_err(ParseError::InvalidInfo)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    parse_genotypes(header, s, record.genotypes_mut())
        .map_err(ParseError::InvalidGenotypes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

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

fn parse_alternate_bases(s: &str) -> io::Result<AlternateBases> {
    match s {
        MISSING => Ok(AlternateBases::default()),
        _ => s
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}
