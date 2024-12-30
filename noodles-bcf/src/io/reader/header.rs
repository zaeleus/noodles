mod format_version;
pub(crate) mod magic_number;
mod vcf_header;

use std::io::{self, BufRead, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{self as vcf, header::StringMaps};

pub(super) use self::{format_version::read_format_version, magic_number::read_magic_number};

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<vcf::Header>
where
    R: Read,
{
    let mut parser = vcf::header::Parser::default();
    let mut string_maps = StringMaps::default();

    let header_len = reader.read_u32::<LittleEndian>().map(u64::from)?;
    let mut header_reader = vcf_header::Reader::new(reader, header_len);

    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf)? != 0 {
        let entry = parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        string_maps
            .insert_entry(&entry)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    header_reader.discard_to_end()?;

    let mut header = parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *header.string_maps_mut() = string_maps;

    Ok(header)
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

        let mut data = 61u32.to_le_bytes().to_vec(); // l_text
        data.extend_from_slice(raw_header);
        data.push(NUL);

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = vcf::Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
