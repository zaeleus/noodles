mod genotypes;
pub mod info;

pub use self::{genotypes::read_genotypes, info::read_info};

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{AlternateBases, Ids, Position, QualityScore, ReferenceBases},
};

use super::value::read_value;
use crate::{
    header::StringMaps,
    lazy::{
        self,
        record::{ChromosomeId, Filters, Value},
    },
};

pub(super) fn read_record<R>(
    reader: &mut R,
    header: &vcf::Header,
    string_maps: &StringMaps,
    buf: &mut Vec<u8>,
    record: &mut vcf::Record,
) -> io::Result<usize>
where
    R: Read,
{
    let l_shared = match reader.read_u32::<LittleEndian>() {
        Ok(n) => usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    let l_indiv = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    buf.resize(l_shared, 0);
    reader.read_exact(buf)?;
    let mut src = &buf[..];
    let (n_fmt, n_sample) = read_site(&mut src, header, string_maps, record)?;

    buf.resize(l_indiv, 0);
    reader.read_exact(buf)?;
    let mut src = &buf[..];

    *record.genotypes_mut() = read_genotypes(
        &mut src,
        header.formats(),
        string_maps.strings(),
        n_sample,
        n_fmt,
    )?;

    Ok(l_shared + l_indiv)
}

fn read_site<R>(
    reader: &mut R,
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &mut vcf::Record,
) -> io::Result<(usize, usize)>
where
    R: Read,
{
    let chrom = read_chrom(reader)?;

    *record.chromosome_mut() = string_maps
        .contigs()
        .get_index(chrom)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid chrom"))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    *record.position_mut() = read_pos(reader)?;

    // TODO
    read_rlen(reader)?;

    *record.quality_score_mut() = read_qual(reader)?;

    let n_info = reader.read_u16::<LittleEndian>().map(usize::from)?;
    let n_allele = reader.read_u16::<LittleEndian>().map(usize::from)?;

    let n_fmt_sample = reader.read_u32::<LittleEndian>()?;
    let n_fmt = usize::from((n_fmt_sample >> 24) as u8);
    let n_sample = usize::try_from(n_fmt_sample & 0xffffff)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *record.ids_mut() = read_id(reader)?;

    let (r#ref, alt) = read_ref_alt(reader, n_allele)?;
    *record.reference_bases_mut() = r#ref;
    *record.alternate_bases_mut() = alt;

    let mut filters = lazy::record::Filters::default();
    read_filter(reader, &mut filters)?;
    *record.filters_mut() = filters.try_into_vcf_record_filters(string_maps.strings())?;

    *record.info_mut() = read_info(reader, header.infos(), string_maps.strings(), n_info)?;

    Ok((n_fmt, n_sample))
}

pub fn read_chrom<R>(reader: &mut R) -> io::Result<ChromosomeId>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>().and_then(|n| {
        ChromosomeId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

pub fn read_rlen<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    reader
        .read_i32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

pub fn read_pos<R>(reader: &mut R) -> io::Result<Position>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n + 1)
            .map(Position::from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

pub fn read_qual<R>(reader: &mut R) -> io::Result<Option<QualityScore>>
where
    R: Read,
{
    use crate::lazy::record::value::Float;

    match reader.read_f32::<LittleEndian>().map(Float::from)? {
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

pub fn read_id<R>(reader: &mut R) -> io::Result<Ids>
where
    R: Read,
{
    match read_value(reader)? {
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

pub fn read_ref_alt<R>(reader: &mut R, len: usize) -> io::Result<(ReferenceBases, AlternateBases)>
where
    R: Read,
{
    let mut alleles = Vec::with_capacity(len);

    for _ in 0..len {
        match read_value(reader)? {
            Some(Value::String(Some(s))) => alleles.push(s),
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

pub fn read_filter<R>(reader: &mut R, filters: &mut Filters) -> io::Result<()>
where
    R: Read,
{
    use super::string_map::read_string_map_indices;

    let filter = filters.as_mut();
    filter.clear();

    let indices = read_string_map_indices(reader)?;
    filter.extend_from_slice(&indices);

    Ok(())
}
