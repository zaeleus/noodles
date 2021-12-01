mod genotypes;
mod site;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

use crate::header::StringMap;

pub fn write_vcf_record<W>(
    writer: &mut W,
    header: &vcf::Header,
    string_map: &StringMap,
    record: &vcf::Record,
) -> io::Result<()>
where
    W: Write,
{
    use self::{genotypes::write_genotypes, site::write_site};

    let mut site_buf = Vec::new();
    write_site(&mut site_buf, header, string_map, record)?;

    let l_shared = u32::try_from(site_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut genotypes_buf = Vec::new();
    let genotypes = record.genotypes();

    if !genotypes.is_empty() {
        write_genotypes(&mut genotypes_buf, string_map, genotypes.keys(), genotypes)?;
    };

    let l_indiv = u32::try_from(genotypes_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32::<LittleEndian>(l_shared)?;
    writer.write_u32::<LittleEndian>(l_indiv)?;
    writer.write_all(&site_buf)?;
    writer.write_all(&genotypes_buf)?;

    Ok(())
}
