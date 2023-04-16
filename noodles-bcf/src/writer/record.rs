mod genotypes;
pub(crate) mod site;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

use crate::header::StringMaps;

pub fn write_record<W>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &vcf::Record,
) -> io::Result<()>
where
    W: Write,
{
    use self::{genotypes::write_genotypes, site::write_site};

    let mut site_buf = Vec::new();
    write_site(&mut site_buf, header, string_maps, record)?;

    let l_shared = u32::try_from(site_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut genotypes_buf = Vec::new();
    let genotypes = record.genotypes();

    if !genotypes.is_empty() {
        write_genotypes(&mut genotypes_buf, header, string_maps.strings(), genotypes)?;
    };

    let l_indiv = u32::try_from(genotypes_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32::<LittleEndian>(l_shared)?;
    writer.write_u32::<LittleEndian>(l_indiv)?;
    writer.write_all(&site_buf)?;
    writer.write_all(&genotypes_buf)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::header::record::value::{map::Contig, Map};

        let header = vcf::Header::builder()
            .add_contig("sq0".parse()?, Map::<Contig>::new())
            .build();

        let string_maps = StringMaps::try_from(&header)?;

        let record = vcf::Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(vcf::record::Position::from(1))
            .set_reference_bases("A".parse()?)
            .build()?;

        let mut buf = Vec::new();
        write_record(&mut buf, &header, &string_maps, &record)?;

        let expected = [
            0x1c, 0x00, 0x00, 0x00, // l_shared = 28
            0x00, 0x00, 0x00, 0x00, // l_indiv = 0
            0x00, 0x00, 0x00, 0x00, // chrom = 0,
            0x00, 0x00, 0x00, 0x00, // pos = 0 (0-based)
            0x01, 0x00, 0x00, 0x00, // rlen = 1
            0x01, 0x00, 0x80, 0x7f, // qual = Float::Missing
            0x00, 0x00, // n_info = 0
            0x01, 0x00, // n_allele = 1
            0x00, // n_fmt = 0
            0x00, 0x00, 0x00, // n_sample = 0
            0x07, // id = None
            0x17, b'A', // ref = A, alt = []
            0x00, // filter = []
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
