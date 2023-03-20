mod chromosome;
mod position;

use std::io;

use self::{chromosome::parse_chromosome, position::parse_position};
use crate::{
    record::{AlternateBases, Filters, Genotypes, Ids, Info, QualityScore, ReferenceBases},
    Header, Record,
};

const MISSING: &str = ".";

pub(super) fn parse_record(mut s: &str, header: &Header, record: &mut Record) -> io::Result<()> {
    let field = next_field(&mut s);
    parse_chromosome(field, record.chromosome_mut())?;

    let field = next_field(&mut s);
    *record.position_mut() = parse_position(field)?;

    let field = next_field(&mut s);
    *record.ids_mut() = parse_ids(field)?;

    let field = next_field(&mut s);
    *record.reference_bases_mut() = parse_reference_bases(field)?;

    let field = next_field(&mut s);
    *record.alternate_bases_mut() = parse_alternate_bases(field)?;

    let field = next_field(&mut s);
    *record.quality_score_mut() = parse_quality_score(field)?;

    let field = next_field(&mut s);
    *record.filters_mut() = parse_filters(field)?;

    let field = next_field(&mut s);
    *record.info_mut() = parse_info(header, field)?;

    *record.genotypes_mut() = parse_genotypes(header, s)?;

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

fn parse_ids(s: &str) -> io::Result<Ids> {
    match s {
        MISSING => Ok(Ids::default()),
        _ => s
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn parse_reference_bases(s: &str) -> io::Result<ReferenceBases> {
    match s {
        MISSING => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing reference bases",
        )),
        _ => s
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
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

fn parse_filters(s: &str) -> io::Result<Option<Filters>> {
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

fn parse_genotypes(header: &Header, s: &str) -> io::Result<Genotypes> {
    if s.is_empty() {
        Ok(Genotypes::default())
    } else {
        Genotypes::parse(s, header).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
