mod genotypes;
pub mod info;

pub use self::{genotypes::read_genotypes, info::read_info};

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::record::{AlternateBases, Ids, Position, QualityScore, ReferenceBases};

use super::value::read_value;
use crate::{
    lazy,
    lazy::record::{ChromosomeId, Filters, Value},
};

pub fn read_lazy_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut lazy::Record,
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

    buf.resize(l_shared, Default::default());
    reader.read_exact(buf)?;
    let mut buf_reader = &buf[..];
    let (n_fmt, n_sample) = read_site(&mut buf_reader, record)?;

    let genotypes = record.genotypes_mut().as_mut();
    genotypes.resize(l_indiv, Default::default());
    reader.read_exact(genotypes)?;
    record.genotypes_mut().set_format_count(n_fmt);
    record.genotypes_mut().set_sample_count(n_sample);

    Ok(l_shared + l_indiv)
}

pub(crate) fn read_site<R>(reader: &mut R, record: &mut lazy::Record) -> io::Result<(usize, usize)>
where
    R: Read,
{
    *record.chromosome_id_mut() = read_chrom(reader)?;
    *record.position_mut() = read_pos(reader)?;

    *record.rlen_mut() = read_rlen(reader)?;

    *record.quality_score_mut() = read_qual(reader)?;

    let n_info = reader.read_u16::<LittleEndian>().map(usize::from)?;
    let n_allele = reader.read_u16::<LittleEndian>().map(usize::from)?;

    let n_fmt_sample = reader.read_u32::<LittleEndian>()?;
    let n_fmt = usize::from((n_fmt_sample >> 24) as u8);
    let n_sample = usize::try_from(n_fmt_sample & 0xffffff)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *record.ids_mut() = read_id(reader)?;

    let (r#ref, alt) = read_ref_alt(reader, n_allele)?;
    record.r#ref = r#ref;
    record.alt = alt;

    read_filter(reader, record.filters_mut())?;

    let info = record.info_mut().as_mut();
    info.clear();
    reader.read_to_end(info)?;
    record.info_mut().set_field_count(n_info);

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

fn read_id<R>(reader: &mut R) -> io::Result<Ids>
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

fn read_ref_alt<R>(reader: &mut R, len: usize) -> io::Result<(ReferenceBases, AlternateBases)>
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

fn read_filter<R>(reader: &mut R, filters: &mut Filters) -> io::Result<()>
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

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub static RAW_HEADER: &str = r#"##fileformat=VCFv4.3
##INFO=<ID=HM3,Number=0,Type=Flag,Description="HM3 membership",IDX=1>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed",IDX=2>
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles called genotypes",IDX=3>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele",IDX=4>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=5>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality",IDX=6>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=7>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele",IDX=8>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer",IDX=9>
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0	sample1	sample2
"#;

    // § Putting it all together (2021-07-27)
    //
    // Note that the data in the reference table mixes big and little endian. INFO string map
    // indices are offset by -79. FORMAT string map indices are offset by +4.
    pub static DATA: [u8; 101] = [
        0x33, 0x00, 0x00, 0x00, // l_shared = 51
        0x2a, 0x00, 0x00, 0x00, // l_indiv = 42
        //
        0x01, 0x00, 0x00, 0x00, // chrom = 1
        0x64, 0x00, 0x00, 0x00, // pos = 100 (base 0)
        0x01, 0x00, 0x00, 0x00, // rlen = 1
        0xcd, 0xcc, 0xf0, 0x41, // qual = 30.1
        0x04, 0x00, // n_info = 4
        0x02, 0x00, // n_allele = 2
        0x03, 0x00, 0x00, // n_sample = 3
        0x05, // n_fmt = 5
        0x57, 0x72, 0x73, 0x31, 0x32, 0x33, // id = "rs123"
        0x17, 0x41, // ref = A
        0x17, 0x43, // alt = C
        0x11, 0x00, // filter = 0 (PASS)
        //
        0x11, 0x01, 0x00, // infos[HM3] = (1, true)
        0x11, 0x02, 0x11, 0x03, // infos[AC] = (2, 3)
        0x11, 0x03, 0x11, 0x06, // infos[AN] = (3, 6)
        0x11, 0x04, 0x17, 0x43, // infos[AA] = (4, "C")
        //
        0x11, 0x05, // formats[GT]
        0x21, // i8[2]
        0x02, 0x02, // 0/0
        0x02, 0x04, // 0/1
        0x04, 0x04, // 1/1
        //
        0x11, 0x06, // formats[GQ]
        0x11, // i8[1]
        0x0a, 0x0a, 0x0a, // [10, 10, 10]
        //
        0x11, 0x07, // formats[DP]
        0x11, // i8[1]
        0x20, 0x30, 0x40, // [32, 48, 64]
        //
        0x11, 0x08, // formats[AD]
        0x21, // i8[2]
        0x20, 0x00, // [32, 0]
        0x20, 0x10, // [32, 16]
        0x00, 0x40, // [0, 64]
        //
        0x11, 0x09, // formats[PL]
        0x31, // i8[3]
        0x00, 0x0a, 0x64, // [0, 10, 100]
        0x0a, 0x00, 0x64, // [10, 0, 100]
        0x64, 0x0a, 0x00, // [100, 10, 0]
    ];

    #[test]
    fn test_read_lazy_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_vcf::{
            header::{format, info},
            record::{
                genotypes::{sample::Value as GenotypeFieldValue, Keys},
                info::field::Value as InfoFieldValue,
                Filters as VcfFilters, Genotypes as VcfGenotypes, Ids, Position,
            },
        };

        use crate::header::StringMaps;

        let header = RAW_HEADER.parse()?;
        let string_maps: StringMaps = RAW_HEADER.parse()?;

        let mut reader = &DATA[..];
        let mut buf = Vec::new();
        let mut record = lazy::Record::default();
        read_lazy_record(&mut reader, &mut buf, &mut record)?;

        assert_eq!(record.chromosome_id(), 1);
        assert_eq!(record.position(), Position::try_from(101)?);
        assert_eq!(record.rlen(), 1);
        assert_eq!(
            record.quality_score(),
            QualityScore::try_from(30.1).map(Some)?
        );
        assert_eq!(record.ids(), &"rs123".parse::<Ids>()?);
        assert_eq!(record.reference_bases(), &"A".parse()?);
        assert_eq!(record.alternate_bases(), &"C".parse()?);

        assert_eq!(
            record
                .filters()
                .try_into_vcf_record_filters(string_maps.strings())?,
            Some(VcfFilters::Pass),
        );

        // info

        let actual = record
            .info()
            .try_into_vcf_record_info(&header, string_maps.strings())?;

        let expected = [
            ("HM3".parse()?, Some(InfoFieldValue::Flag)),
            (info::key::ALLELE_COUNT, Some(InfoFieldValue::Integer(3))),
            (
                info::key::TOTAL_ALLELE_COUNT,
                Some(InfoFieldValue::Integer(6)),
            ),
            (
                info::key::ANCESTRAL_ALLELE,
                Some(InfoFieldValue::String(String::from("C"))),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        // genotypes

        let actual = record
            .genotypes()
            .try_into_vcf_record_genotypes(&header, string_maps.strings())?;

        let expected = VcfGenotypes::new(
            Keys::try_from(vec![
                format::key::GENOTYPE,
                format::key::CONDITIONAL_GENOTYPE_QUALITY,
                format::key::READ_DEPTH,
                format::key::READ_DEPTHS,
                format::key::ROUNDED_GENOTYPE_LIKELIHOODS,
            ])?,
            vec![
                vec![
                    Some(GenotypeFieldValue::String(String::from("0/0"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(32)),
                    Some(GenotypeFieldValue::IntegerArray(vec![Some(32), Some(0)])),
                    Some(GenotypeFieldValue::IntegerArray(vec![
                        Some(0),
                        Some(10),
                        Some(100),
                    ])),
                ],
                vec![
                    Some(GenotypeFieldValue::String(String::from("0/1"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(48)),
                    Some(GenotypeFieldValue::IntegerArray(vec![Some(32), Some(16)])),
                    Some(GenotypeFieldValue::IntegerArray(vec![
                        Some(10),
                        Some(0),
                        Some(100),
                    ])),
                ],
                vec![
                    Some(GenotypeFieldValue::String(String::from("1/1"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(64)),
                    Some(GenotypeFieldValue::IntegerArray(vec![Some(0), Some(64)])),
                    Some(GenotypeFieldValue::IntegerArray(vec![
                        Some(100),
                        Some(10),
                        Some(0),
                    ])),
                ],
            ],
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
