mod genotypes;
mod site;

pub use self::{
    genotypes::read_genotypes,
    site::{read_site, Site},
};

use std::io::{self, Read};

use noodles_vcf::{self as vcf, record::genotypes::Genotype};

use crate::header::StringMap;

pub fn read_record<R>(
    reader: &mut R,
    header: &vcf::Header,
    string_map: &StringMap,
) -> io::Result<(Site, Vec<Genotype>)>
where
    R: Read,
{
    let site = read_site(reader, header, string_map)?;

    let genotypes = if site.n_sample == 0 {
        Vec::new()
    } else {
        read_genotypes(
            reader,
            string_map,
            site.n_sample as usize,
            usize::from(site.n_fmt),
        )?
    };

    Ok((site, genotypes))
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use noodles_vcf::record::{Filters, Info};

    use crate::record::value::Float;

    use super::*;

    #[test]
    fn test_read_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_vcf::record;

        // § Putting it all together (2021-01-13)
        //
        // Note that the data in the reference table mixes big and little endian. INFO string map
        // indices are offset by -79. FORMAT string map indices are offset by +4.
        let data = [
            0x01, 0x00, 0x00, 0x00, // chrom = 1
            0x64, 0x00, 0x00, 0x00, // pos = 100 (base 0)
            0x01, 0x00, 0x00, 0x00, // rlen = 1
            0xcd, 0xcc, 0xf0, 0x41, // qual = 30.1
            0x04, 0x00, 0x02, 0x00, // n_allele_info (allele count, info count) = (2, 4)
            0x03, 0x00, 0x00, 0x05, // n_fmt_sample (format count, sample_count) = (5, 3)
            0x57, 0x72, 0x73, 0x31, 0x32, 0x33, // id = "rs123"
            0x17, 0x41, // ref = A
            0x17, 0x43, // alt = C
            0x11, 0x00, // filter = 0 (PASS)
            //
            0x11, 0x01, 0x11, 0x01, // infos[HM3] = (1, true)
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

        let mut reader = &data[..];

        let raw_header = r#"##fileformat=VCFv4.3
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

        let header = raw_header.parse()?;
        let string_map = raw_header.parse()?;

        let (actual_site, actual_genotypes) = read_record(&mut reader, &header, &string_map)?;

        let expected_site = Site {
            chrom: 1,
            pos: 100,
            rlen: 1,
            qual: Float::from(30.1),
            n_info: 4,
            n_allele: 2,
            n_sample: 3,
            n_fmt: 5,
            id: "rs123".parse()?,
            ref_alt: vec![String::from("A"), String::from("C")],
            filter: Filters::try_from_iter(&["PASS"])?,
            info: Info::try_from(vec![
                record::info::Field::new("HM3".parse()?, record::info::field::Value::Flag),
                record::info::Field::new(
                    record::info::field::Key::AlleleCount,
                    record::info::field::Value::Integer(3),
                ),
                record::info::Field::new(
                    record::info::field::Key::TotalAlleleCount,
                    record::info::field::Value::Integer(6),
                ),
                record::info::Field::new(
                    record::info::field::Key::AncestralAllele,
                    record::info::field::Value::String(String::from("C")),
                ),
            ])?,
        };

        assert_eq!(actual_site, expected_site);

        let expected_genotypes = vec![
            Genotype::try_from(vec![
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::Genotype,
                    Some(record::genotypes::genotype::field::Value::String(
                        String::from("0/0"),
                    )),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ConditionalGenotypeQuality,
                    Some(record::genotypes::genotype::field::Value::Integer(10)),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ReadDepth,
                    Some(record::genotypes::genotype::field::Value::Integer(32)),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ReadDepths,
                    Some(record::genotypes::genotype::field::Value::IntegerArray(
                        vec![Some(32), Some(0)],
                    )),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::RoundedGenotypeLikelihoods,
                    Some(record::genotypes::genotype::field::Value::IntegerArray(
                        vec![Some(0), Some(10), Some(100)],
                    )),
                ),
            ])?,
            Genotype::try_from(vec![
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::Genotype,
                    Some(record::genotypes::genotype::field::Value::String(
                        String::from("0/1"),
                    )),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ConditionalGenotypeQuality,
                    Some(record::genotypes::genotype::field::Value::Integer(10)),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ReadDepth,
                    Some(record::genotypes::genotype::field::Value::Integer(48)),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ReadDepths,
                    Some(record::genotypes::genotype::field::Value::IntegerArray(
                        vec![Some(32), Some(16)],
                    )),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::RoundedGenotypeLikelihoods,
                    Some(record::genotypes::genotype::field::Value::IntegerArray(
                        vec![Some(10), Some(0), Some(100)],
                    )),
                ),
            ])?,
            Genotype::try_from(vec![
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::Genotype,
                    Some(record::genotypes::genotype::field::Value::String(
                        String::from("1/1"),
                    )),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ConditionalGenotypeQuality,
                    Some(record::genotypes::genotype::field::Value::Integer(10)),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ReadDepth,
                    Some(record::genotypes::genotype::field::Value::Integer(64)),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::ReadDepths,
                    Some(record::genotypes::genotype::field::Value::IntegerArray(
                        vec![Some(0), Some(64)],
                    )),
                ),
                record::genotypes::genotype::Field::new(
                    record::genotypes::genotype::field::Key::RoundedGenotypeLikelihoods,
                    Some(record::genotypes::genotype::field::Value::IntegerArray(
                        vec![Some(100), Some(10), Some(0)],
                    )),
                ),
            ])?,
        ];

        assert_eq!(actual_genotypes, expected_genotypes);

        Ok(())
    }
}
