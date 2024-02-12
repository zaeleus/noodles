use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::Record;

pub fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
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

    let site_buf = record.fields_mut().site_buf_mut();
    site_buf.resize(l_shared, 0);
    reader.read_exact(site_buf)?;

    record.fields_mut().index()?;

    let samples_buf = record.fields_mut().samples_buf_mut();
    samples_buf.resize(l_indiv, 0);
    reader.read_exact(samples_buf)?;

    Ok(l_shared + l_indiv)
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
    fn test_read_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_core::Position;
        use noodles_vcf::variant::record_buf::{
            info::{self, field::Value as InfoFieldValue},
            samples::{
                self,
                sample::{value::Array, Value as GenotypeFieldValue},
                Keys,
            },
            Samples as VcfGenotypes,
        };

        use crate::header::StringMaps;

        let header = RAW_HEADER.parse()?;
        let string_maps: StringMaps = RAW_HEADER.parse()?;

        let mut reader = &DATA[..];
        let mut record = Record::default();
        read_record(&mut reader, &mut record)?;

        assert_eq!(record.chromosome_id()?, 1);
        assert_eq!(
            record.position().transpose()?,
            Some(Position::try_from(101)?)
        );
        assert_eq!(record.rlen()?, 1);
        assert_eq!(record.quality_score()?, Some(30.1));
        assert_eq!(record.ids().as_ref(), b"rs123");
        assert_eq!(record.reference_bases().as_ref(), b"A");
        assert_eq!(
            record
                .alternate_bases()
                .iter()
                .collect::<Result<Vec<_>, _>>()?,
            [Some("C")]
        );

        assert_eq!(
            record
                .filters()
                .iter(&string_maps)?
                .collect::<io::Result<Vec<_>>>()?,
            ["PASS"],
        );

        // info

        let actual = record
            .info()
            .try_into_vcf_record_info(&header, &string_maps)?;

        let expected = [
            (String::from("HM3"), Some(InfoFieldValue::Flag)),
            (
                String::from(info::field::key::ALLELE_COUNT),
                Some(InfoFieldValue::from(3)),
            ),
            (
                String::from(info::field::key::TOTAL_ALLELE_COUNT),
                Some(InfoFieldValue::from(6)),
            ),
            (
                String::from(info::field::key::ANCESTRAL_ALLELE),
                Some(InfoFieldValue::from("C")),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        // genotypes

        let actual = record
            .samples()?
            .try_into_vcf_record_samples(&header, string_maps.strings())?;

        let expected = VcfGenotypes::new(
            Keys::try_from(vec![
                String::from(samples::keys::key::GENOTYPE),
                String::from(samples::keys::key::CONDITIONAL_GENOTYPE_QUALITY),
                String::from(samples::keys::key::READ_DEPTH),
                String::from(samples::keys::key::READ_DEPTHS),
                String::from(samples::keys::key::ROUNDED_GENOTYPE_LIKELIHOODS),
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
