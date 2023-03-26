mod chromosome;
mod filters;
mod genotypes;
mod ids;
mod position;
mod reference_bases;

use std::io;

use self::{
    chromosome::parse_chromosome, filters::parse_filters, genotypes::parse_genotypes,
    ids::parse_ids, position::parse_position, reference_bases::parse_reference_bases,
};
use crate::{
    record::{AlternateBases, Info, QualityScore},
    Header, Record,
};

const MISSING: &str = ".";

pub(super) fn parse_record(mut s: &str, header: &Header, record: &mut Record) -> io::Result<()> {
    let field = next_field(&mut s);
    parse_chromosome(field, record.chromosome_mut())?;

    let field = next_field(&mut s);
    *record.position_mut() = parse_position(field)?;

    let field = next_field(&mut s);
    parse_ids(field, record.ids_mut())?;

    let field = next_field(&mut s);
    parse_reference_bases(field, record.reference_bases_mut())?;

    let field = next_field(&mut s);
    *record.alternate_bases_mut() = parse_alternate_bases(field)?;

    let field = next_field(&mut s);
    *record.quality_score_mut() = parse_quality_score(field)?;

    let field = next_field(&mut s);
    parse_filters(field, record.filters_mut())?;

    let field = next_field(&mut s);
    *record.info_mut() = parse_info(header, field)?;

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

fn parse_quality_score(s: &str) -> io::Result<Option<QualityScore>> {
    match s {
        MISSING => Ok(None),
        _ => s
            .parse()
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn parse_info(header: &Header, s: &str) -> io::Result<Info> {
    match s {
        MISSING => Ok(Info::default()),
        _ => Info::try_from_str(s, header.infos())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}
