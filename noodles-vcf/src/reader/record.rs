mod chromosome;
mod filters;
mod genotypes;
mod ids;
mod info;
mod position;
mod reference_bases;

use std::io;

use self::{
    chromosome::parse_chromosome, filters::parse_filters, genotypes::parse_genotypes,
    ids::parse_ids, info::parse_info, position::parse_position,
    reference_bases::parse_reference_bases,
};
use crate::{
    record::{AlternateBases, QualityScore},
    Header, Record,
};

const MISSING: &str = ".";

pub(super) fn parse_record(mut s: &str, header: &Header, record: &mut Record) -> io::Result<()> {
    let field = next_field(&mut s);
    parse_chromosome(field, record.chromosome_mut())?;

    let field = next_field(&mut s);
    *record.position_mut() = parse_position(field)?;

    record.ids_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_ids(field, record.ids_mut())?;
    }

    let field = next_field(&mut s);
    parse_reference_bases(field, record.reference_bases_mut())?;

    let field = next_field(&mut s);
    *record.alternate_bases_mut() = parse_alternate_bases(field)?;

    let field = next_field(&mut s);
    *record.quality_score_mut() = match field {
        MISSING => None,
        _ => parse_quality_score(field).map(Some)?,
    };

    let field = next_field(&mut s);
    match field {
        MISSING => {
            record.filters_mut().take();
        }
        _ => parse_filters(field, record.filters_mut())?,
    }

    record.info_mut().clear();
    let field = next_field(&mut s);
    if field != MISSING {
        parse_info(header, field, record.info_mut())?;
    }

    parse_genotypes(header, s, record.genotypes_mut())?;

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

fn parse_quality_score(s: &str) -> io::Result<QualityScore> {
    s.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
