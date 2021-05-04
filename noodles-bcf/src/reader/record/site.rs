use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{Filters, Ids, Info},
};

use crate::{
    header::StringMap,
    reader::value::{read_value, Float, Value},
};

#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub chrom: i32,
    pub pos: i32,
    pub rlen: i32,
    pub qual: Float,
    pub n_info: u16,
    pub n_allele: u16,
    pub n_sample: u32,
    pub n_fmt: u8,
    pub id: Ids,
    pub ref_alt: Vec<String>,
    pub filter: Filters,
    pub info: Info,
}

pub fn read_site<R>(
    reader: &mut R,
    header: &vcf::Header,
    string_map: &StringMap,
) -> io::Result<Site>
where
    R: Read,
{
    let chrom = reader.read_i32::<LittleEndian>()?;
    let pos = reader.read_i32::<LittleEndian>()?;

    let rlen = reader.read_i32::<LittleEndian>()?;

    let qual = reader.read_f32::<LittleEndian>().map(Float::from)?;

    let n_info = reader.read_u16::<LittleEndian>()?;
    let n_allele = reader.read_u16::<LittleEndian>()?;

    let n_fmt_sample = reader.read_u32::<LittleEndian>()?;
    let n_fmt = (n_fmt_sample >> 24) as u8;
    let n_sample = n_fmt_sample & 0xffffff;

    let id = read_id(reader)?;
    let ref_alt = read_ref_alt(reader, usize::from(n_allele))?;
    let filter = read_filter(reader, string_map)?;
    let info = read_info(reader, header.infos(), string_map, usize::from(n_info))?;

    Ok(Site {
        chrom,
        pos,
        rlen,
        qual,
        n_info,
        n_allele,
        n_sample,
        n_fmt,
        id,
        ref_alt,
        filter,
        info,
    })
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
