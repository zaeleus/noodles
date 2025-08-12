mod bases;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_sequence_id;

use std::io::{self, Write};

use noodles_vcf::{
    self as vcf,
    header::StringMaps,
    variant::{Record, record::AlternateBases},
};

use self::{
    bases::write_bases, filters::write_filters, ids::write_ids, info::write_info,
    position::write_position, quality_score::write_quality_score,
    reference_sequence_id::write_reference_sequence_id,
};
use crate::io::writer::num::{write_i32_le, write_u16_le, write_u32_le};

const MAX_SAMPLE_NAME_COUNT: u32 = (1 << 24) - 1;

pub fn write_site<W, R>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &R,
) -> io::Result<()>
where
    W: Write,
    R: Record + ?Sized,
{
    let reference_sequence_name = record.reference_sequence_name(header)?;
    write_reference_sequence_id(writer, string_maps.contigs(), reference_sequence_name)?;

    let position = record.variant_start().transpose()?;
    write_position(writer, position)?;

    let span = record.variant_span(header)?;
    write_rlen(writer, span)?;

    let quality_score = record.quality_score().transpose()?;
    write_quality_score(writer, quality_score)?;

    let n_info = u16::try_from(record.info().as_ref().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u16_le(writer, n_info)?;

    write_n_allele(writer, record.alternate_bases().len())?;

    let samples = record.samples()?;
    write_n_fmt_sample(
        writer,
        header.sample_names().len(),
        samples.column_names(header).count(),
    )?;

    write_ids(writer, record.ids())?;
    write_bases(writer, record.reference_bases(), record.alternate_bases())?;
    write_filters(writer, header, string_maps, record.filters())?;
    write_info(writer, header, string_maps, record.info())?;

    Ok(())
}

pub(crate) fn write_rlen<W>(writer: &mut W, span: usize) -> io::Result<()>
where
    W: Write,
{
    let rlen = i32::try_from(span).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, rlen)
}

pub(crate) fn write_n_allele<W>(writer: &mut W, alternate_base_count: usize) -> io::Result<()>
where
    W: Write,
{
    const REFERENCE_BASE_COUNT: usize = 1;

    let n_allele = alternate_base_count
        .checked_add(REFERENCE_BASE_COUNT)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "attempt to add with overflow"))
        .and_then(|n| {
            u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    write_u16_le(writer, n_allele)?;

    Ok(())
}

pub(crate) fn write_n_fmt_sample<W>(
    writer: &mut W,
    sample_count: usize,
    format_count: usize,
) -> io::Result<()>
where
    W: Write,
{
    let n_sample =
        u32::try_from(sample_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if n_sample > MAX_SAMPLE_NAME_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid sample name count: {n_sample}"),
        ));
    }

    let n_fmt =
        u8::try_from(format_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let n_fmt_sample = (u32::from(n_fmt) << 24) | n_sample;
    write_u32_le(writer, n_fmt_sample)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_rlen() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        write_rlen(&mut buf, 6)?;
        assert_eq!(buf, &[0x06, 0x00, 0x00, 0x00]);
        Ok(())
    }

    #[test]
    fn test_n_allele() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, alternate_base_count: usize, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_n_allele(buf, alternate_base_count)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x01, 0x00])?;
        t(&mut buf, 1, &[0x02, 0x00])?;
        t(&mut buf, usize::from(u16::MAX - 1), &[0xff, 0xff])?;

        buf.clear();
        assert!(matches!(
            write_n_allele(&mut buf, usize::from(u16::MAX)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        assert!(matches!(
            write_n_allele(&mut buf, usize::MAX),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_n_fmt_sample() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            sample_count: usize,
            format_count: usize,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_n_fmt_sample(buf, sample_count, format_count)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, 0, &[0x00, 0x00, 0x00, 0x00])?;
        t(&mut buf, 1, 0, &[0x01, 0x00, 0x00, 0x00])?;
        t(&mut buf, 0, 1, &[0x00, 0x00, 0x00, 0x01])?;

        t(&mut buf, 8, 13, &[0x08, 0x00, 0x00, 0x0d])?;

        buf.clear();
        assert!(matches!(
            write_n_fmt_sample(&mut buf, 1 << 24, 0),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        assert!(matches!(
            write_n_fmt_sample(&mut buf, 0, 1 << 8),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
