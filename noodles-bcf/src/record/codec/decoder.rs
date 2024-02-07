mod bases;
mod chromosome_id;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod raw_value;
mod samples;
mod string_map;
mod value;

use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf as vcf;

use self::info::read_info;
pub(crate) use self::{
    bases::read_ref_alt, chromosome_id::read_chrom, filters::read_filter, ids::read_id,
    position::read_pos, quality_score::read_qual, string_map::read_string_map_entry,
};
pub use self::{samples::read_samples, value::read_value};
use crate::header::StringMaps;

pub fn read_site(
    src: &mut &[u8],
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &mut vcf::variant::RecordBuf,
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

    let filter_ids = read_filter(src)?;

    let raw_filters: Vec<_> = filter_ids
        .into_iter()
        .map(|i| {
            string_maps.strings().get_index(i).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid string map index: {i}"),
                )
            })
        })
        .collect::<io::Result<_>>()?;

    *record.filters_mut() = if raw_filters.is_empty() {
        None
    } else {
        vcf::variant::record_buf::Filters::try_from_iter(raw_filters)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
    };

    read_info(
        src,
        header.infos(),
        string_maps.strings(),
        n_info,
        record.info_mut(),
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((n_fmt, n_sample))
}

pub fn read_rlen(src: &mut &[u8]) -> io::Result<usize> {
    src.read_i32::<LittleEndian>()
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}
