mod info;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

use super::value::write_value;
use crate::{
    header::{
        string_maps::{ContigStringMap, StringStringMap},
        StringMaps,
    },
    lazy::record::value::{Float, Value},
};

use self::info::write_info;

const MAX_SAMPLE_NAME_COUNT: u32 = (1 << 24) - 1;

pub fn write_site<W>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &vcf::Record,
) -> io::Result<()>
where
    W: Write,
{
    write_chrom(writer, string_maps.contigs(), record.chromosome())?;
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
    write_filter(writer, string_maps.strings(), record.filters())?;
    write_info(writer, string_maps.strings(), record.info())?;

    Ok(())
}

fn write_chrom<W>(
    writer: &mut W,
    contig_string_map: &ContigStringMap,
    chromosome: &str,
) -> io::Result<()>
where
    W: Write,
{
    let chrom = contig_string_map
        .get_index_of(chromosome)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("chromosome not in string map: {chromosome}"),
            )
        })
        .and_then(|i| {
            i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    writer.write_i32::<LittleEndian>(chrom)
}

pub(crate) fn write_pos<W>(writer: &mut W, position: vcf::record::Position) -> io::Result<()>
where
    W: Write,
{
    let pos = i32::try_from(usize::from(position))
        .map(|n| n - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_i32::<LittleEndian>(pos)
}

pub(crate) fn write_rlen<W>(
    writer: &mut W,
    start: vcf::record::Position,
    end: vcf::record::Position,
) -> io::Result<()>
where
    W: Write,
{
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

pub(crate) fn write_qual<W>(
    writer: &mut W,
    quality_score: Option<vcf::record::QualityScore>,
) -> io::Result<()>
where
    W: Write,
{
    let float = quality_score
        .map(|qs| Float::from(f32::from(qs)))
        .unwrap_or(Float::Missing);

    writer.write_f32::<LittleEndian>(f32::from(float))
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

pub(crate) fn write_id<W>(writer: &mut W, ids: &vcf::record::Ids) -> io::Result<()>
where
    W: Write,
{
    if ids.is_empty() {
        let value = Some(Value::String(None));
        write_value(writer, value)
    } else {
        let s = ids.to_string();
        let value = Some(Value::String(Some(&s)));
        write_value(writer, value)
    }
}

pub(crate) fn write_ref_alt<W>(
    writer: &mut W,
    reference_bases: &vcf::record::ReferenceBases,
    alternate_bases: &vcf::record::AlternateBases,
) -> io::Result<()>
where
    W: Write,
{
    let r#ref = reference_bases.to_string();
    let ref_value = Some(Value::String(Some(&r#ref)));
    write_value(writer, ref_value)?;

    if !alternate_bases.is_empty() {
        for allele in alternate_bases.iter() {
            let alt = allele.to_string();
            let alt_value = Some(Value::String(Some(&alt)));
            write_value(writer, alt_value)?;
        }
    }

    Ok(())
}

fn write_filter<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    filters: Option<&vcf::record::Filters>,
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::Filters;

    use crate::record::codec::encoder::string_map::write_string_map_indices;

    let indices = match filters {
        None => Vec::new(),
        Some(Filters::Pass) => vec![0],
        Some(Filters::Fail(ids)) => ids
            .iter()
            .map(|id| {
                string_string_map.get_index_of(id).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("filter missing from string map: {id}"),
                    )
                })
            })
            .collect::<Result<_, _>>()?,
    };

    write_string_map_indices(writer, &indices)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_chrom() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            contig_string_map: &ContigStringMap,
            chromosome: &str,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_chrom(buf, contig_string_map, chromosome)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let mut contig_string_map = ContigStringMap::default();
        contig_string_map.insert("sq0".into());
        contig_string_map.insert("sq1".into());

        t(
            &mut buf,
            &contig_string_map,
            "sq0",
            &[0x00, 0x00, 0x00, 0x00],
        )?;

        t(
            &mut buf,
            &contig_string_map,
            "sq1",
            &[0x01, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_pos() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::record::Position;

        fn t(buf: &mut Vec<u8>, position: Position, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_pos(buf, position)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Position::from(0), &[0xff, 0xff, 0xff, 0xff])?;
        t(&mut buf, Position::from(1), &[0x00, 0x00, 0x00, 0x00])?;
        t(
            &mut buf,
            Position::from((1 << 31) - 1),
            &[0xfe, 0xff, 0xff, 0x7f],
        )?;

        buf.clear();
        assert!(matches!(
            write_pos(&mut buf, Position::from(1 << 32)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_rlen() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::record::Position;

        fn t(buf: &mut Vec<u8>, start: Position, end: Position, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_rlen(buf, start, end)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            Position::from(8),
            Position::from(13),
            &[0x06, 0x00, 0x00, 0x00],
        )?;

        buf.clear();
        assert!(matches!(
            write_rlen(&mut buf, Position::from(13), Position::from(8)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_write_qual() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::record::QualityScore;

        fn t(
            buf: &mut Vec<u8>,
            quality_score: Option<vcf::record::QualityScore>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_qual(buf, quality_score)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[0x01, 0x00, 0x80, 0x7f])?;

        t(
            &mut buf,
            QualityScore::try_from(0.0).map(Some)?,
            &[0x00, 0x00, 0x00, 0x00],
        )?;

        t(
            &mut buf,
            QualityScore::try_from(8.0).map(Some)?,
            &[0x00, 0x00, 0x00, 0x41],
        )?;

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

    #[test]
    fn test_write_id() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::record::Ids;

        fn t(buf: &mut Vec<u8>, ids: &Ids, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_id(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Ids::default(), &[0x07])?;
        t(&mut buf, &"nd0".parse()?, &[0x37, b'n', b'd', b'0'])?;
        t(
            &mut buf,
            &"nd0;nd1".parse()?,
            &[0x77, b'n', b'd', b'0', b';', b'n', b'd', b'1'],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_ref_alt() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::record::{AlternateBases, ReferenceBases};

        fn t(
            buf: &mut Vec<u8>,
            reference_bases: &ReferenceBases,
            alternate_bases: &AlternateBases,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_ref_alt(buf, reference_bases, alternate_bases)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &"A".parse()?,
            &AlternateBases::default(),
            &[0x17, b'A'],
        )?;

        t(
            &mut buf,
            &"A".parse()?,
            &"G".parse()?,
            &[0x17, b'A', 0x17, b'G'],
        )?;

        t(
            &mut buf,
            &"A".parse()?,
            &"G,T".parse()?,
            &[0x17, b'A', 0x17, b'G', 0x17, b'T'],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_filter() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::{
            header::record::value::{map::Filter, Map},
            record::Filters,
        };

        fn t(
            buf: &mut Vec<u8>,
            string_map: &StringStringMap,
            filters: Option<&Filters>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_filter(buf, string_map, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = vcf::Header::builder()
            .add_filter("PASS", Map::<Filter>::pass())
            .add_filter(
                "s50",
                Map::<Filter>::new("Less than 50% of samples have data"),
            )
            .add_filter("q10", Map::<Filter>::new("Quality below 10"))
            .build();

        let string_maps = StringMaps::try_from(&header)?;

        let mut buf = Vec::new();

        t(&mut buf, string_maps.strings(), None, &[0x00])?;

        let filters = Filters::Pass;
        t(
            &mut buf,
            string_maps.strings(),
            Some(&filters),
            &[0x11, 0x00],
        )?;

        let filters = Filters::try_from_iter(["q10"])?;
        t(
            &mut buf,
            string_maps.strings(),
            Some(&filters),
            &[0x11, 0x02],
        )?;

        let filters = Filters::try_from_iter(["q10", "s50"])?;
        t(
            &mut buf,
            string_maps.strings(),
            Some(&filters),
            &[0x21, 0x02, 0x01],
        )?;

        Ok(())
    }
}
