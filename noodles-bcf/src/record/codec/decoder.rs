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

use std::{io, num::NonZero};

use noodles_vcf as vcf;

use self::info::read_info;
pub(crate) use self::{
    bases::read_ref_alt, chromosome_id::read_chrom, filters::read_filter, ids::read_id,
    position::read_pos, quality_score::read_qual, string_map::read_string_map_entry,
};
pub use self::{samples::read_samples, value::read_value};
use crate::io::reader::num::{read_i32_le, read_u16_le, read_u32_le};

pub fn read_site(
    src: &mut &[u8],
    header: &vcf::Header,
    record: &mut vcf::variant::RecordBuf,
) -> io::Result<(usize, usize)> {
    let chrom = read_chrom(src)?;

    *record.reference_sequence_name_mut() = header
        .string_maps()
        .contigs()
        .get_index(chrom)
        .map(String::from)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid chrom"))?;

    *record.variant_start_mut() = read_pos(src)?;

    // TODO
    read_rlen(src)?;

    *record.quality_score_mut() = read_qual(src)?;

    let n_info = read_info_field_count(src)?;
    let n_allele = read_allele_count(src)?;
    let (n_fmt, n_sample) = read_key_and_sample_counts(src)?;

    *record.ids_mut() = read_id(src)?;

    let (r#ref, alt) = read_ref_alt(src, n_allele)?;
    *record.reference_bases_mut() = r#ref;
    *record.alternate_bases_mut() = alt;

    let filter_ids = read_filter(src)?;

    *record.filters_mut() = filter_ids
        .into_iter()
        .map(|i| {
            header
                .string_maps()
                .strings()
                .get_index(i)
                .map(String::from)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid string map index: {i}"),
                    )
                })
        })
        .collect::<io::Result<_>>()?;

    read_info(src, header, n_info, record.info_mut())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((n_fmt, n_sample))
}

fn read_rlen(src: &mut &[u8]) -> io::Result<usize> {
    read_i32_le(src)
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

fn read_info_field_count(src: &mut &[u8]) -> io::Result<usize> {
    read_u16_le(src).map(usize::from)
}

fn read_allele_count(src: &mut &[u8]) -> io::Result<NonZero<usize>> {
    read_u16_le(src).and_then(|n| {
        NonZero::try_from(usize::from(n)).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_key_and_sample_counts(src: &mut &[u8]) -> io::Result<(usize, usize)> {
    let n_fmt_sample = read_u32_le(src)?;

    let n_fmt = usize::from((n_fmt_sample >> 24) as u8);

    let n_sample = usize::try_from(n_fmt_sample & 0xffffff)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((n_fmt, n_sample))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_allele_count() -> io::Result<()> {
        assert_eq!(
            read_allele_count(&mut &[0x01, 0x00][..])?,
            NonZero::<usize>::MIN
        );

        assert_eq!(
            read_allele_count(&mut &[0xff, 0xff][..])?,
            const { NonZero::new(u16::MAX as usize).unwrap() }
        );

        assert!(matches!(
            read_allele_count(&mut &[0x00, 0x00][..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
