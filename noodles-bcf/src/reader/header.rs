use std::io::{self, BufRead, BufReader, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf as vcf;

use crate::header::StringMaps;

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<(vcf::Header, StringMaps)>
where
    R: Read,
{
    const NUL: u8 = 0x00;

    let l_text = reader.read_u32::<LittleEndian>().map(u64::from)?;

    let mut parser = vcf::header::Parser::default();
    let mut string_maps = StringMaps::default();

    let mut header_reader = BufReader::new(reader.take(l_text));
    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf)? != 0 {
        if buf == [NUL] {
            break;
        }

        let entry = parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        string_maps
            .insert_entry(&entry)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    let header = parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((header, string_maps))
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() -> io::Result<()> {
        use vcf::header::FileFormat;

        const NUL: u8 = 0x00;

        let raw_header = b"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
";

        let mut data = 60u32.to_le_bytes().to_vec(); // l_text = 22
        data.extend_from_slice(raw_header);
        data.push(NUL);

        let mut reader = &data[..];
        let (actual_header, actual_string_maps) = read_header(&mut reader)?;

        let expected_header = vcf::Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .build();

        let expected_string_maps = StringMaps::default();

        assert_eq!(actual_header, expected_header);
        assert_eq!(actual_string_maps, expected_string_maps);

        Ok(())
    }
}
