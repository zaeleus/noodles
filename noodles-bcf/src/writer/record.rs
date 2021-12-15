use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::{record::Filters, Record};

pub(super) fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    let mut shared = Vec::new();
    write_site(&mut shared, record)?;
    let l_shared =
        u32::try_from(shared.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let indiv = record.genotypes().as_ref();
    let l_indiv =
        u32::try_from(indiv.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32::<LittleEndian>(l_shared)?;
    writer.write_u32::<LittleEndian>(l_indiv)?;
    writer.write_all(&shared)?;
    writer.write_all(indiv)?;

    Ok(())
}

fn write_site<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    use super::vcf_record::site::{
        write_id, write_n_allele, write_n_fmt_sample, write_pos, write_qual, write_ref_alt,
        write_rlen,
    };

    let chrom = record.chromosome_id();
    writer.write_i32::<LittleEndian>(chrom)?;

    write_pos(writer, record.position())?;

    let end = record.end()?;
    write_rlen(writer, record.position(), end)?;

    write_qual(writer, record.quality_score())?;

    let n_info = u16::try_from(record.info().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_info)?;

    write_n_allele(writer, record.alternate_bases().len())?;

    let genotypes = record.genotypes();
    write_n_fmt_sample(writer, genotypes.len(), genotypes.format_count())?;

    write_id(writer, record.ids())?;
    write_ref_alt(writer, record.reference_bases(), record.alternate_bases())?;

    write_filter(writer, record.filters())?;

    writer.write_all(record.info().as_ref())?;

    Ok(())
}

fn write_filter<W>(writer: &mut W, filters: &Filters) -> io::Result<()>
where
    W: Write,
{
    use crate::writer::value::write_value;

    if filters.is_empty() {
        write_value(writer, None)?;
    } else {
        todo!();
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut buf = Vec::new();
        let record = Record::default();
        write_record(&mut buf, &record)?;

        let expected = [
            0x1d, 0x00, 0x00, 0x00, // l_shared = 29
            0x00, 0x00, 0x00, 0x00, // l_indiv = 0
            0x00, 0x00, 0x00, 0x00, // chrom = 0,
            0x00, 0x00, 0x00, 0x00, // pos = 0 (0-based)
            0x01, 0x00, 0x00, 0x00, // rlen = 1
            0x01, 0x00, 0x80, 0x7f, // qual = Float::Missing
            0x00, 0x00, // n_info = 0
            0x01, 0x00, // n_allele = 1
            0x00, 0x00, 0x00, // n_sample = 0
            0x00, // n_fmt = 0
            0x07, // id = None
            0x17, b'A', // ref = [A]
            0x07, // alt = []
            0x00, // filter = []
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
