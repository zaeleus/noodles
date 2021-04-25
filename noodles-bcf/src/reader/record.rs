use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{Filters, Ids, Info, QualityScore},
};

use crate::header::StringMap;

use super::value::{read_value, Value};

#[allow(dead_code, unused_variables)]
pub fn read_site<R>(reader: &mut R, header: &vcf::Header, string_map: &StringMap) -> io::Result<()>
where
    R: Read,
{
    let chrom = reader.read_i32::<LittleEndian>()?;
    let pos = reader.read_i32::<LittleEndian>()?;

    let rlen = reader.read_i32::<LittleEndian>()?;

    let qual = reader.read_f32::<LittleEndian>().and_then(|value| {
        QualityScore::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let n_allele_info = reader.read_i32::<LittleEndian>()?;
    let allele_count = (n_allele_info >> 16) as u16;
    let info_count = (n_allele_info & 0xffff) as u16;

    let n_fmt_sample = reader.read_i32::<LittleEndian>()?;
    let format_count = (n_fmt_sample >> 24) as u8;
    let sample_count = n_fmt_sample & 0xffffff;

    let id = read_id(reader)?;
    let ref_alt = read_ref_alt(reader, usize::from(allele_count))?;
    let filter = read_filter(reader, string_map)?;
    let info = read_info(reader, header.infos(), string_map, usize::from(info_count))?;

    Ok(())
}

fn read_id<R>(reader: &mut R) -> io::Result<Ids>
where
    R: Read,
{
    match read_value(reader) {
        Ok(Value::String(Some(id))) => id
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        Ok(Value::String(None)) => Ok(Ids::default()),
        Ok(v) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected string, got {:?}", v),
        )),
        Err(e) => Err(e),
    }
}

fn read_ref_alt<R>(reader: &mut R, len: usize) -> io::Result<Vec<String>>
where
    R: Read,
{
    let mut alleles = Vec::with_capacity(len);

    for _ in 0..len {
        match read_value(reader) {
            Ok(Value::String(Some(s))) => alleles.push(s),
            Ok(Value::String(None)) => alleles.push(String::from(".")),
            Ok(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected string, got {:?}", v),
                ))
            }
            Err(e) => return Err(e),
        }
    }

    Ok(alleles)
}

fn read_filter<R>(reader: &mut R, string_map: &StringMap) -> io::Result<Filters>
where
    R: Read,
{
    let indices = match read_value(reader) {
        Ok(Value::Int8(None)) => Vec::new(),
        Ok(Value::Int8(Some(i))) => vec![i],
        Ok(Value::Int8Array(indices)) => indices,
        Ok(v) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ))
        }
        Err(e) => return Err(e),
    };

    let raw_filters: Vec<_> = indices
        .iter()
        .map(|&i| {
            usize::try_from(i)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .and_then(|j| {
                    string_map.get(j).ok_or_else(|| {
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
    use vcf::{
        header::info::ty::Type,
        record::info::{self, Field},
    };

    let mut fields = Vec::with_capacity(len);

    for _ in 0..len {
        let key: info::field::Key = match read_value(reader)? {
            Value::Int8(Some(i)) => usize::try_from(i)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .and_then(|j| {
                    string_map.get(j).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("invalid string map index: {}", j),
                        )
                    })
                })
                .and_then(|s| {
                    s.parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })?,
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected i8, got {:?}", v),
                ));
            }
        };

        let info = infos.get(&key).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing header INFO record for {}", key),
            )
        })?;

        let value = match info.ty() {
            Type::Integer => match read_value(reader)? {
                Value::Int8(Some(n)) => info::field::Value::Integer(i32::from(n)),
                Value::Int8Array(values) => {
                    info::field::Value::IntegerArray(values.into_iter().map(i32::from).collect())
                }
                Value::Int16(Some(n)) => info::field::Value::Integer(i32::from(n)),
                Value::Int16Array(values) => {
                    info::field::Value::IntegerArray(values.into_iter().map(i32::from).collect())
                }
                Value::Int32(Some(n)) => info::field::Value::Integer(n),
                Value::Int32Array(values) => info::field::Value::IntegerArray(values),
                v => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("type mismatch: expected {}, got {:?}", Type::Integer, v),
                    ));
                }
            },
            Type::Flag => match read_value(reader)? {
                Value::Int8(Some(1)) => info::field::Value::Flag,
                v => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("type mismatch: expected {}, got {:?}", Type::Flag, v),
                    ));
                }
            },
            Type::String => match read_value(reader)? {
                Value::String(Some(s)) => info::field::Value::String(s),
                v => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("type mismatch: expected {}, got {:?}", Type::String, v),
                    ))
                }
            },
            ty => todo!("unhandled INFO value type: {:?}", ty),
        };

        fields.push(Field::new(key, value));
    }

    Info::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_site() -> Result<(), Box<dyn std::error::Error>> {
        // § Putting it all together (2021-01-13)
        //
        // Note that the data in the reference table mixes big and little endian. INFO string map
        // indices are offset by -79.
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
            0x11, 0x01, 0x11, 0x01, // infos[HM3] = (1, true)
            0x11, 0x02, 0x11, 0x03, // infos[AC] = (2, 3)
            0x11, 0x03, 0x11, 0x06, // infos[AN] = (3, 6)
            0x11, 0x04, 0x17, 0x43, // infos[AA] = (4, "C")
        ];

        let mut reader = &data[..];

        let raw_header = r#"##fileformat=VCFv4.3
##INFO=<ID=HM3,Number=0,Type=Flag,Description="HM3 membership",IDX=1>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed",IDX=2>
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles called genotypes",IDX=3>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele",IDX=4>
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let header = raw_header.parse()?;
        let string_map = raw_header.parse()?;

        read_site(&mut reader, &header, &string_map)?;

        Ok(())
    }
}
