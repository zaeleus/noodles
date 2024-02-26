use std::io::{self, Write};

use noodles_vcf::header::string_maps::ContigStringMap;

use byteorder::{LittleEndian, WriteBytesExt};

pub(super) fn write_reference_sequence_id<W>(
    writer: &mut W,
    contig_string_map: &ContigStringMap,
    chromosome: &str,
) -> io::Result<()>
where
    W: Write,
{
    let chrom = contig_string_map
        .get_index_of(chromosome)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("chromosome not in string map: {chromosome}"),
            )
        })
        .and_then(|i| {
            i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })?;

    writer.write_i32::<LittleEndian>(chrom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_chrom() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            contig_string_map: &ContigStringMap,
            chromosome: &str,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_reference_sequence_id(buf, contig_string_map, chromosome)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let mut contig_string_map = ContigStringMap::default();
        contig_string_map.insert("sq0".into());
        contig_string_map.insert("sq1".into());

        t(
            &mut buf,
            &contig_string_map,
            "sq0",
            &[0x00, 0x00, 0x00, 0x00],
        )?;

        t(
            &mut buf,
            &contig_string_map,
            "sq1",
            &[0x01, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }
}
