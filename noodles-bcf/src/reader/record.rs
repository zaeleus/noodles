mod site;

pub use self::site::Site;

use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{Filters, Genotype, Ids, Info},
};

use crate::header::StringMap;

use super::value::{read_type, read_value, Float, Value};

#[allow(dead_code)]
pub fn read_site<R>(
    reader: &mut R,
    header: &vcf::Header,
    string_map: &StringMap,
) -> io::Result<(Site, Vec<Genotype>)>
where
    R: Read,
{
    let chrom = reader.read_i32::<LittleEndian>()?;
    let pos = reader.read_i32::<LittleEndian>()?;

    let rlen = reader.read_i32::<LittleEndian>()?;

    let qual = reader.read_f32::<LittleEndian>().map(Float::from)?;

    let n_allele_info = reader.read_i32::<LittleEndian>()?;
    let allele_count = (n_allele_info >> 16) as u16;
    let info_count = (n_allele_info & 0xffff) as u16;

    let n_fmt_sample = reader.read_u32::<LittleEndian>()?;
    let format_count = (n_fmt_sample >> 24) as u8;
    let sample_count = n_fmt_sample & 0xffffff;

    let id = read_id(reader)?;
    let ref_alt = read_ref_alt(reader, usize::from(allele_count))?;
    let filter = read_filter(reader, string_map)?;
    let info = read_info(reader, header.infos(), string_map, usize::from(info_count))?;

    let site = Site {
        chrom,
        pos,
        rlen,
        qual,
        n_allele_info,
        n_fmt_sample,
        id,
        ref_alt,
        filter,
        info,
    };

    let genotypes = read_genotypes(
        reader,
        string_map,
        sample_count as usize,
        usize::from(format_count),
    )?;

    Ok((site, genotypes))
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
            format!("expected string, got {:?}", v),
        )),
    }
}

fn read_ref_alt<R>(reader: &mut R, len: usize) -> io::Result<Vec<String>>
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
                    format!("expected string, got {:?}", v),
                ))
            }
        }
    }

    Ok(alleles)
}

fn read_filter<R>(reader: &mut R, string_map: &StringMap) -> io::Result<Filters>
where
    R: Read,
{
    use crate::reader::value::Int8;

    let indices = match read_value(reader)? {
        Some(Value::Int8(None)) | None => Vec::new(),
        Some(Value::Int8(Some(Int8::Value(i)))) => vec![i],
        Some(Value::Int8Array(indices)) => indices,
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ))
        }
    };

    let raw_filters: Vec<_> = indices
        .iter()
        .map(|&i| {
            usize::try_from(i)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .and_then(|j| {
                    string_map.get_index(j).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("invalid string map index: {}", j),
                        )
                    })
                })
        })
        .collect::<Result<_, _>>()?;

    Filters::try_from_iter(raw_filters).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_info<R>(
    reader: &mut R,
    infos: &vcf::header::Infos,
    string_map: &StringMap,
    len: usize,
) -> io::Result<Info>
where
    R: Read,
{
    use vcf::record::info::Field;

    let mut fields = Vec::with_capacity(len);

    for _ in 0..len {
        let key = read_info_key(reader, string_map)?;

        let info = infos.get(&key).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing header INFO record for {}", key),
            )
        })?;

        let value = read_info_value(reader, &info)?;

        let field = Field::new(key, value);
        fields.push(field);
    }

    Info::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_info_key<R>(
    reader: &mut R,
    string_map: &StringMap,
) -> io::Result<vcf::record::info::field::Key>
where
    R: Read,
{
    use crate::reader::value::Int8;

    match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(i)))) => usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|j| {
                string_map.get_index(j).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid string map index: {}", j),
                    )
                })
            })
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ));
        }
    }
}

fn read_info_value<R>(
    reader: &mut R,
    info: &vcf::header::Info,
) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    use vcf::{header::info::Type, record::info};

    use crate::reader::value::{Int16, Int32, Int8};

    let value = match info.ty() {
        Type::Integer => match read_value(reader)? {
            Some(Value::Int8(Some(Int8::Value(n)))) => info::field::Value::Integer(i32::from(n)),
            Some(Value::Int8Array(values)) => {
                info::field::Value::IntegerArray(values.into_iter().map(i32::from).collect())
            }
            Some(Value::Int16(Some(Int16::Value(n)))) => info::field::Value::Integer(i32::from(n)),
            Some(Value::Int16Array(values)) => {
                info::field::Value::IntegerArray(values.into_iter().map(i32::from).collect())
            }
            Some(Value::Int32(Some(Int32::Value(n)))) => info::field::Value::Integer(n),
            Some(Value::Int32Array(values)) => info::field::Value::IntegerArray(values),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Integer, v),
                ));
            }
        },
        Type::Flag => match read_value(reader)? {
            Some(Value::Int8(Some(Int8::Value(1)))) | None => info::field::Value::Flag,
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Flag, v),
                ));
            }
        },
        Type::Float => match read_value(reader)? {
            Some(Value::Float(Some(Float::Value(n)))) => info::field::Value::Float(n),
            Some(Value::FloatArray(values)) => info::field::Value::FloatArray(values),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Float, v),
                ))
            }
        },
        Type::String => match read_value(reader)? {
            Some(Value::String(Some(s))) => info::field::Value::String(s),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::String, v),
                ))
            }
        },
        ty => todo!("unhandled INFO value type: {:?}", ty),
    };

    Ok(value)
}

fn read_genotypes<R>(
    reader: &mut R,
    string_map: &StringMap,
    sample_count: usize,
    format_count: usize,
) -> io::Result<Vec<Genotype>>
where
    R: Read,
{
    use vcf::record::genotype::{self, Field};

    let mut genotypes = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_genotype_key(reader, string_map)?;

        let values = if key == genotype::field::Key::Genotype {
            read_genotype_genotype_values(reader, sample_count)?
        } else {
            read_genotype_values(reader, sample_count)?
        };

        for (fields, value) in genotypes.iter_mut().zip(values) {
            let field = Field::new(key.clone(), value);
            fields.push(field);
        }
    }

    genotypes
        .into_iter()
        .map(|fields| {
            Genotype::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .collect()
}

fn read_genotype_key<R>(
    reader: &mut R,
    string_map: &StringMap,
) -> io::Result<vcf::record::genotype::field::Key>
where
    R: Read,
{
    use crate::reader::value::Int8;

    match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(i)))) => usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|j| {
                string_map.get_index(j).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid string map index: {}", j),
                    )
                })
            })
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ));
        }
    }
}

fn read_genotype_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    use super::value::{Int8, Type};

    let mut values = Vec::with_capacity(sample_count);

    match read_type(reader)? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
            1 => {
                for _ in 0..sample_count {
                    let value = reader.read_i8().map(Int8::from)?;

                    match value {
                        Int8::Value(n) => {
                            values.push(Some(genotype::field::Value::Integer(i32::from(n))))
                        }
                        Int8::Missing => values.push(None),
                        _ => todo!("unhandled i8 value: {:?}", value),
                    }
                }
            }
            _ => {
                for _ in 0..sample_count {
                    let mut buf = vec![0; len];
                    reader.read_i8_into(&mut buf)?;
                    let value = genotype::field::Value::IntegerArray(
                        buf.into_iter()
                            .map(Int8::from)
                            .map(|value| match value {
                                Int8::Value(n) => Some(i32::from(n)),
                                Int8::Missing => None,
                                _ => todo!("unhanlded i8 array value: {:?}", value),
                            })
                            .collect(),
                    );
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn read_genotype_genotype_values<R>(
    reader: &mut R,
    sample_count: usize,
) -> io::Result<Vec<Option<vcf::record::genotype::field::Value>>>
where
    R: Read,
{
    use vcf::record::genotype;

    use super::value::Type;

    let mut values = Vec::with_capacity(sample_count);

    match read_type(reader)? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
            1 => {
                for _ in 0..sample_count {
                    let value = reader
                        .read_i8()
                        .map(|v| parse_genotype_genotype_values(&[v]))
                        .map(genotype::field::Value::String)?;

                    values.push(Some(value));
                }
            }
            _ => {
                for _ in 0..sample_count {
                    let mut buf = vec![0; len];
                    reader.read_i8_into(&mut buf)?;
                    let value =
                        genotype::field::Value::String(parse_genotype_genotype_values(&buf));
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_genotype_values(values: &[i8]) -> String {
    use crate::reader::value::Int8;

    let mut genotype = String::new();

    for (i, &value) in values.iter().enumerate() {
        if let Int8::EndOfVector = Int8::from(value) {
            break;
        }

        let j = (value >> 1) - 1;
        let is_phased = value & 0x01 == 1;

        if i > 0 {
            if is_phased {
                genotype.push('|');
            } else {
                genotype.push('/');
            }
        }

        if j == -1 {
            genotype.push('.');
        } else {
            genotype.push_str(&format!("{}", j));
        }
    }

    genotype
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_site() -> Result<(), Box<dyn std::error::Error>> {
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

        let (actual_site, actual_genotypes) = read_site(&mut reader, &header, &string_map)?;

        let expected_site = Site {
            chrom: 1,
            pos: 100,
            rlen: 1,
            qual: Float::from(30.1),
            n_allele_info: 2 << 16 | 4,
            n_fmt_sample: 5 << 24 | 3,
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
                record::genotype::Field::new(
                    record::genotype::field::Key::Genotype,
                    Some(record::genotype::field::Value::String(String::from("0/0"))),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ConditionalGenotypeQuality,
                    Some(record::genotype::field::Value::Integer(10)),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ReadDepth,
                    Some(record::genotype::field::Value::Integer(32)),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ReadDepths,
                    Some(record::genotype::field::Value::IntegerArray(vec![
                        Some(32),
                        Some(0),
                    ])),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::RoundedGenotypeLikelihoods,
                    Some(record::genotype::field::Value::IntegerArray(vec![
                        Some(0),
                        Some(10),
                        Some(100),
                    ])),
                ),
            ])?,
            Genotype::try_from(vec![
                record::genotype::Field::new(
                    record::genotype::field::Key::Genotype,
                    Some(record::genotype::field::Value::String(String::from("0/1"))),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ConditionalGenotypeQuality,
                    Some(record::genotype::field::Value::Integer(10)),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ReadDepth,
                    Some(record::genotype::field::Value::Integer(48)),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ReadDepths,
                    Some(record::genotype::field::Value::IntegerArray(vec![
                        Some(32),
                        Some(16),
                    ])),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::RoundedGenotypeLikelihoods,
                    Some(record::genotype::field::Value::IntegerArray(vec![
                        Some(10),
                        Some(0),
                        Some(100),
                    ])),
                ),
            ])?,
            Genotype::try_from(vec![
                record::genotype::Field::new(
                    record::genotype::field::Key::Genotype,
                    Some(record::genotype::field::Value::String(String::from("1/1"))),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ConditionalGenotypeQuality,
                    Some(record::genotype::field::Value::Integer(10)),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ReadDepth,
                    Some(record::genotype::field::Value::Integer(64)),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::ReadDepths,
                    Some(record::genotype::field::Value::IntegerArray(vec![
                        Some(0),
                        Some(64),
                    ])),
                ),
                record::genotype::Field::new(
                    record::genotype::field::Key::RoundedGenotypeLikelihoods,
                    Some(record::genotype::field::Value::IntegerArray(vec![
                        Some(100),
                        Some(10),
                        Some(0),
                    ])),
                ),
            ])?,
        ];

        assert_eq!(actual_genotypes, expected_genotypes);

        Ok(())
    }

    #[test]
    fn test_parse_genotype_genotype_values() {
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x02]), "0/0");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04]), "0/1");
        assert_eq!(parse_genotype_genotype_values(&[0x04, 0x04]), "1/1");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x05]), "0|1");
        assert_eq!(parse_genotype_genotype_values(&[0x00, 0x00]), "./.");
        assert_eq!(parse_genotype_genotype_values(&[0x02]), "0");
        assert_eq!(parse_genotype_genotype_values(&[0x04]), "1");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04, 0x06]), "0/1/2");
        assert_eq!(parse_genotype_genotype_values(&[0x02, 0x04, 0x07]), "0/1|2");
        assert_eq!(parse_genotype_genotype_values(&[0x02, -127]), "0");
    }
}
