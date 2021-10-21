mod info;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

use crate::{
    header::StringMap,
    record::value::{Float, Value},
    writer::value::write_value,
};

use self::info::write_info;

const MAX_SAMPLE_NAME_COUNT: u32 = (1 << 24) - 1;

pub fn write_site<W>(
    writer: &mut W,
    header: &vcf::Header,
    string_map: &StringMap,
    record: &vcf::Record,
) -> io::Result<()>
where
    W: Write,
{
    write_chrom(writer, header.contigs(), record.chromosome())?;
    write_pos(writer, record.position())?;

    let start = i32::from(record.position());
    let end = record
        .end()
        .map(i32::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let rlen = end - start + 1;
    writer.write_i32::<LittleEndian>(rlen)?;

    write_qual(writer, record.quality_score())?;

    let n_info = u16::try_from(record.info().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_info)?;

    let alternate_bases_len = if record.alternate_bases().is_empty() {
        1
    } else {
        record.alternate_bases().len()
    };

    let n_allele = u16::try_from(1 + alternate_bases_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_allele)?;

    let n_sample = u32::try_from(header.sample_names().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if n_sample > MAX_SAMPLE_NAME_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid sample name count: {}", n_sample),
        ));
    }

    let n_fmt = record
        .format()
        .map(|format| {
            u8::try_from(format.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })
        .unwrap_or(Ok(0))?;

    let n_fmt_sample = u32::from(n_fmt) << 24 | n_sample;
    writer.write_u32::<LittleEndian>(n_fmt_sample)?;

    write_id(writer, record.ids())?;
    write_ref_alt(writer, record.reference_bases(), record.alternate_bases())?;
    write_filter(writer, string_map, record.filters())?;
    write_info(writer, string_map, record.info())?;

    Ok(())
}

fn write_chrom<W>(
    writer: &mut W,
    contigs: &vcf::header::Contigs,
    chromosome: &vcf::record::Chromosome,
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::Chromosome;

    let chrom = match chromosome {
        Chromosome::Name(name) => contigs
            .get_index_of(name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("chromosome not in string map: {}", name),
                )
            })
            .and_then(|i| {
                i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?,
        Chromosome::Symbol(_) => todo!("unhandled chromosome: {:?}", chromosome),
    };

    writer.write_i32::<LittleEndian>(chrom)
}

fn write_pos<W>(writer: &mut W, position: vcf::record::Position) -> io::Result<()>
where
    W: Write,
{
    let pos = i32::from(position) - 1;
    writer.write_i32::<LittleEndian>(pos)
}

fn write_qual<W>(writer: &mut W, quality_score: vcf::record::QualityScore) -> io::Result<()>
where
    W: Write,
{
    let float = quality_score.map(Float::from).unwrap_or(Float::Missing);
    writer.write_f32::<LittleEndian>(f32::from(float))
}

fn write_id<W>(writer: &mut W, ids: &vcf::record::Ids) -> io::Result<()>
where
    W: Write,
{
    let value = if ids.is_empty() {
        Some(Value::String(None))
    } else {
        Some(Value::String(Some(ids.to_string())))
    };

    write_value(writer, value)
}

fn write_ref_alt<W>(
    writer: &mut W,
    reference_bases: &vcf::record::ReferenceBases,
    alternate_bases: &vcf::record::AlternateBases,
) -> io::Result<()>
where
    W: Write,
{
    let r#ref = reference_bases.to_string();
    let ref_value = Some(Value::String(Some(r#ref)));
    write_value(writer, ref_value)?;

    if alternate_bases.is_empty() {
        write_value(writer, Some(Value::String(None)))?;
    } else {
        for allele in alternate_bases.iter() {
            let alt_value = Some(Value::String(Some(allele.to_string())));
            write_value(writer, alt_value)?;
        }
    }

    Ok(())
}

fn write_filter<W>(
    writer: &mut W,
    string_map: &StringMap,
    filters: &vcf::record::Filters,
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::Filters;

    let indices = match filters {
        Filters::Missing => Vec::new(),
        Filters::Pass => vec![0],
        Filters::Fail(ids) => ids
            .iter()
            .map(|id| {
                string_map
                    .get_index_of(id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("filter missing from string map: {}", id),
                        )
                    })
                    .and_then(|i| {
                        i8::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
                    })
            })
            .collect::<Result<_, _>>()?,
    };

    if indices.is_empty() {
        write_value(writer, None)
    } else {
        let value = Some(Value::Int8Array(indices));
        write_value(writer, value)
    }
}
