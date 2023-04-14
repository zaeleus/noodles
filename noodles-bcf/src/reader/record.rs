mod genotypes;
pub mod info;

pub use self::{genotypes::read_genotypes, info::read_info};

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::record::{AlternateBases, Ids, Position, QualityScore, ReferenceBases};

use super::value::read_value;
use crate::lazy::record::{ChromosomeId, Filters, Value};

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

pub fn read_id<R>(reader: &mut R) -> io::Result<Ids>
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

pub fn read_ref_alt<R>(reader: &mut R, len: usize) -> io::Result<(ReferenceBases, AlternateBases)>
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

pub fn read_filter<R>(reader: &mut R, filters: &mut Filters) -> io::Result<()>
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
