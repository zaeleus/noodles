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

    let end = record
        .end()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_rlen(writer, record.position(), end)?;

    write_qual(writer, record.quality_score())?;

    let n_info = u16::try_from(record.info().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_info)?;

    write_n_allele(writer, record.alternate_bases().len())?;

    write_n_fmt_sample(
        writer,
        header.sample_names().len(),
        record.genotypes().keys().len(),
    )?;

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

fn write_rlen<W>(
    writer: &mut W,
    start: vcf::record::Position,
    end: vcf::record::Position,
) -> io::Result<()>
where
    W: Write,
{
    let rlen = i32::from(start) - i32::from(end) + 1;
    writer.write_i32::<LittleEndian>(rlen)
}

fn write_qual<W>(writer: &mut W, quality_score: Option<vcf::record::QualityScore>) -> io::Result<()>
where
    W: Write,
{
    let float = quality_score
        .map(|qs| Float::from(f32::from(qs)))
        .unwrap_or(Float::Missing);

    writer.write_f32::<LittleEndian>(f32::from(float))
}

fn write_n_allele<W>(writer: &mut W, alternate_base_count: usize) -> io::Result<()>
where
    W: Write,
{
    use std::cmp;

    let n = cmp::max(1, alternate_base_count);
    let n_allele =
        u16::try_from(1 + n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_allele)?;

    Ok(())
}

fn write_n_fmt_sample<W>(writer: &mut W, sample_count: usize, format_count: usize) -> io::Result<()>
where
    W: Write,
{
    let n_sample =
        u32::try_from(sample_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if n_sample > MAX_SAMPLE_NAME_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid sample name count: {}", n_sample),
        ));
    }

    let n_fmt =
        u8::try_from(format_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let n_fmt_sample = u32::from(n_fmt) << 24 | n_sample;
    writer.write_u32::<LittleEndian>(n_fmt_sample)?;

    Ok(())
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
    filters: Option<&vcf::record::Filters>,
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::Filters;

    let indices = match filters {
        None => Vec::new(),
        Some(Filters::Pass) => vec![0],
        Some(Filters::Fail(ids)) => ids
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
