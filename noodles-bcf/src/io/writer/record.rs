use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf::{self as vcf, header::StringMaps};

pub fn write_record<W>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    record: &vcf::variant::RecordBuf,
) -> io::Result<()>
where
    W: Write,
{
    use crate::record::codec::encoder::{samples::write_samples, site::write_site};

    let mut site_buf = Vec::new();
    write_site(&mut site_buf, header, string_maps, record)?;

    let l_shared = u32::try_from(site_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut samples_buf = Vec::new();
    let samples = record.samples();

    if !samples.is_empty() {
        write_samples(&mut samples_buf, header, string_maps.strings(), samples)?;
    };

    let l_indiv = u32::try_from(samples_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32::<LittleEndian>(l_shared)?;
    writer.write_u32::<LittleEndian>(l_indiv)?;
    writer.write_all(&site_buf)?;
    writer.write_all(&samples_buf)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_core::Position;
        use vcf::header::record::value::{map::Contig, Map};

        let header = vcf::Header::builder()
            .add_contig("sq0", Map::<Contig>::new())
            .build();

        let string_maps = StringMaps::try_from(&header)?;

        let record = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_position(Position::MIN)
            .set_reference_bases("A")
            .build();

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
