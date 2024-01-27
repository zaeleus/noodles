use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::lazy;
use crate::record::codec::decoder::{
    read_chrom, read_filter, read_id, read_pos, read_qual, read_ref_alt, read_rlen,
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

    let genotypes = record.genotypes.as_mut();
    genotypes.resize(l_indiv, Default::default());
    reader.read_exact(genotypes)?;
    record.genotypes.set_format_count(n_fmt);
    record.genotypes.set_sample_count(n_sample);

    Ok(l_shared + l_indiv)
}

pub(crate) fn read_site(src: &mut &[u8], record: &mut lazy::Record) -> io::Result<(usize, usize)> {
    record.chrom = read_chrom(src)?;
    record.pos = read_pos(src)?;

    record.rlen = read_rlen(src)?;

    record.qual = read_qual(src)?;

    let n_info = src.read_u16::<LittleEndian>().map(usize::from)?;
    let n_allele = src.read_u16::<LittleEndian>().map(usize::from)?;

    let n_fmt_sample = src.read_u32::<LittleEndian>()?;
    let n_fmt = usize::from((n_fmt_sample >> 24) as u8);
    let n_sample = usize::try_from(n_fmt_sample & 0xffffff)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    record.id = read_id(src)?;

    let (r#ref, alt) = read_ref_alt(src, n_allele)?;
    record.r#ref = r#ref;
    record.alt = alt;

    read_filter(src, &mut record.filter)?;

    let info = record.info.as_mut();
    info.clear();
    src.read_to_end(info)?;
    record.info.set_field_count(n_info);

    Ok((n_fmt, n_sample))
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
        use noodles_vcf::record::{
            genotypes::{
                self,
                sample::{value::Array, Value as GenotypeFieldValue},
                Keys,
            },
            info::{self, field::Value as InfoFieldValue},
            Filters as VcfFilters, Genotypes as VcfGenotypes, Ids, Position,
        };

        use crate::header::StringMaps;

        let header = RAW_HEADER.parse()?;
        let string_maps: StringMaps = RAW_HEADER.parse()?;

        let mut reader = &DATA[..];
        let mut buf = Vec::new();
        let mut record = lazy::Record::default();
        read_lazy_record(&mut reader, &mut buf, &mut record)?;

        assert_eq!(record.chromosome_id(), 1);
        assert_eq!(record.position(), Position::from(101));
        assert_eq!(record.rlen(), 1);
        assert_eq!(record.quality_score(), Some(30.1));
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
            (
                info::field::key::ALLELE_COUNT,
                Some(InfoFieldValue::from(3)),
            ),
            (
                info::field::key::TOTAL_ALLELE_COUNT,
                Some(InfoFieldValue::from(6)),
            ),
            (
                info::field::key::ANCESTRAL_ALLELE,
                Some(InfoFieldValue::from("C")),
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
                genotypes::keys::key::GENOTYPE,
                genotypes::keys::key::CONDITIONAL_GENOTYPE_QUALITY,
                genotypes::keys::key::READ_DEPTH,
                genotypes::keys::key::READ_DEPTHS,
                genotypes::keys::key::ROUNDED_GENOTYPE_LIKELIHOODS,
            ])?,
            vec![
                vec![
                    Some(GenotypeFieldValue::String(String::from("0/0"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(32)),
                    Some(GenotypeFieldValue::Array(Array::Integer(vec![
                        Some(32),
                        Some(0),
                    ]))),
                    Some(GenotypeFieldValue::Array(Array::Integer(vec![
                        Some(0),
                        Some(10),
                        Some(100),
                    ]))),
                ],
                vec![
                    Some(GenotypeFieldValue::String(String::from("0/1"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(48)),
                    Some(GenotypeFieldValue::Array(Array::Integer(vec![
                        Some(32),
                        Some(16),
                    ]))),
                    Some(GenotypeFieldValue::Array(Array::Integer(vec![
                        Some(10),
                        Some(0),
                        Some(100),
                    ]))),
                ],
                vec![
                    Some(GenotypeFieldValue::String(String::from("1/1"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(64)),
                    Some(GenotypeFieldValue::Array(Array::Integer(vec![
                        Some(0),
                        Some(64),
                    ]))),
                    Some(GenotypeFieldValue::Array(Array::Integer(vec![
                        Some(100),
                        Some(10),
                        Some(0),
                    ]))),
                ],
            ],
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
