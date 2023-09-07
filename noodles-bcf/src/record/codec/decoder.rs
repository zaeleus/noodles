mod genotypes;
pub mod info;
mod string_map;
mod value;

pub use self::{genotypes::read_genotypes, info::read_info};

use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{AlternateBases, Ids, Position, QualityScore, ReferenceBases},
};

use self::value::read_value;
use crate::{
    header::StringMaps,
    lazy::{
        self,
        record::{ChromosomeId, Filters, Value},
    },
};

pub fn read_site(
    src: &mut &[u8],
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &mut vcf::Record,
) -> io::Result<(usize, usize)> {
    let chrom = read_chrom(src)?;

    *record.chromosome_mut() = string_maps
        .contigs()
        .get_index(chrom)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid chrom"))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    *record.position_mut() = read_pos(src)?;

    // TODO
    read_rlen(src)?;

    *record.quality_score_mut() = read_qual(src)?;

    let n_info = src.read_u16::<LittleEndian>().map(usize::from)?;
    let n_allele = src.read_u16::<LittleEndian>().map(usize::from)?;

    let n_fmt_sample = src.read_u32::<LittleEndian>()?;
    let n_fmt = usize::from((n_fmt_sample >> 24) as u8);
    let n_sample = usize::try_from(n_fmt_sample & 0xffffff)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *record.ids_mut() = read_id(src)?;

    let (r#ref, alt) = read_ref_alt(src, n_allele)?;
    *record.reference_bases_mut() = r#ref;
    *record.alternate_bases_mut() = alt;

    let mut filters = lazy::record::Filters::default();
    read_filter(src, &mut filters)?;
    *record.filters_mut() = filters.try_into_vcf_record_filters(string_maps.strings())?;

    *record.info_mut() = read_info(src, header.infos(), string_maps.strings(), n_info)?;

    Ok((n_fmt, n_sample))
}

pub fn read_chrom(src: &mut &[u8]) -> io::Result<ChromosomeId> {
    src.read_i32::<LittleEndian>().and_then(|n| {
        ChromosomeId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

pub fn read_rlen(src: &mut &[u8]) -> io::Result<usize> {
    src.read_i32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

pub fn read_pos(src: &mut &[u8]) -> io::Result<Position> {
    src.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n + 1)
            .map(Position::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

pub fn read_qual(src: &mut &[u8]) -> io::Result<Option<QualityScore>> {
    use crate::lazy::record::value::Float;

    match src.read_f32::<LittleEndian>().map(Float::from)? {
        Float::Value(value) => QualityScore::try_from(value)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
        Float::Missing => Ok(None),
        qual => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid qual: {qual:?}"),
        )),
    }
}

pub fn read_id(src: &mut &[u8]) -> io::Result<Ids> {
    match read_value(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
        Some(Value::String(Some(id))) => id
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        Some(Value::String(None)) => Ok(Ids::default()),
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid id: expected string, got {v:?}"),
        )),
    }
}

pub fn read_ref_alt(src: &mut &[u8], len: usize) -> io::Result<(ReferenceBases, AlternateBases)> {
    let mut alleles = Vec::with_capacity(len);

    for _ in 0..len {
        match read_value(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
            Some(Value::String(Some(s))) => alleles.push(s.into()),
            Some(Value::String(None)) => alleles.push(String::from(".")),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid ref_alt: expected string, got {v:?}"),
                ))
            }
        }
    }

    let (raw_reference_bases, raw_alternate_bases) = alleles.split_at(1);

    let reference_bases = raw_reference_bases
        .first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing reference bases"))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    let alternate_bases = raw_alternate_bases
        .iter()
        .map(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(AlternateBases::from)?;

    Ok((reference_bases, alternate_bases))
}

pub fn read_filter(reader: &mut &[u8], filters: &mut Filters) -> io::Result<()> {
    use self::string_map::read_string_map_indices;

    let filter = filters.as_mut();
    filter.clear();

    let indices = read_string_map_indices(reader)?;
    filter.extend_from_slice(&indices);

    Ok(())
}
