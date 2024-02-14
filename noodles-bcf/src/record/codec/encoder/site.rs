mod bases;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_sequence_id;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_core::Position;
use noodles_vcf as vcf;

use self::{
    bases::write_bases, filters::write_filters, ids::write_ids, info::write_info,
    position::write_position, quality_score::write_quality_score,
    reference_sequence_id::write_reference_sequence_id,
};
use crate::header::StringMaps;

const MAX_SAMPLE_NAME_COUNT: u32 = (1 << 24) - 1;

pub fn write_site<W>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &vcf::variant::RecordBuf,
) -> io::Result<()>
where
    W: Write,
{
    write_reference_sequence_id(
        writer,
        string_maps.contigs(),
        record.reference_sequence_name(),
    )?;

    write_position(writer, record.position())?;

    let end = record
        .end()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_rlen(writer, record.position(), end)?;

    write_quality_score(writer, record.quality_score())?;

    let n_info = u16::try_from(record.info().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_info)?;

    write_n_allele(writer, record.alternate_bases().len())?;

    write_n_fmt_sample(
        writer,
        header.sample_names().len(),
        record.samples().keys().len(),
    )?;

    write_ids(writer, record.ids())?;
    write_bases(writer, record.reference_bases(), record.alternate_bases())?;
    write_filters(writer, string_maps.strings(), record.filters())?;
    write_info(writer, string_maps.strings(), record.info())?;

    Ok(())
}

pub(crate) fn write_rlen<W>(
    writer: &mut W,
    start: Option<Position>,
    end: Position,
) -> io::Result<()>
where
    W: Write,
{
    let Some(start) = start else {
        todo!();
    };

    if start > end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("record start position ({start}) > end position ({end})"),
        ));
    }

    let rlen = i32::try_from(usize::from(end) - usize::from(start) + 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_i32::<LittleEndian>(rlen)
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

    writer.write_u16::<LittleEndian>(n_allele)?;

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

    let n_fmt_sample = u32::from(n_fmt) << 24 | n_sample;
    writer.write_u32::<LittleEndian>(n_fmt_sample)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_rlen() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            start: Option<Position>,
            end: Position,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_rlen(buf, start, end)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            Some(Position::try_from(8)?),
            Position::try_from(13)?,
            &[0x06, 0x00, 0x00, 0x00],
        )?;

        buf.clear();
        assert!(matches!(
            write_rlen(&mut buf, Some(Position::try_from(13)?), Position::try_from(8)?),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

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
