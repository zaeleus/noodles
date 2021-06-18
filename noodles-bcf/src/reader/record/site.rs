mod info;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{
    self as vcf,
    record::{Filters, Ids, Info},
};

use crate::{
    header::StringMap,
    reader::{string_map::read_string_map_indices, value::read_value},
    record::{value::Float, Value},
};

use self::info::read_info;

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
    let indices = read_string_map_indices(reader)?;

    let raw_filters: Vec<_> = indices
        .iter()
        .map(|&i| {
            string_map.get_index(i).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid string map index: {}", i),
                )
            })
        })
        .collect::<Result<_, _>>()?;

    Filters::try_from_iter(raw_filters).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
