use std::{
    ffi::CStr,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf as vcf;

use crate::header::StringMaps;

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<(vcf::Header, StringMaps)>
where
    R: Read,
{
    let raw_header = read_raw_header(reader)?;

    let header = raw_header
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let string_maps = StringMaps::from(&header);

    Ok((header, string_maps))
}

pub fn read_raw_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    let l_text = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; l_text];
    reader.read_exact(&mut buf)?;

    CStr::from_bytes_with_nul(&buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|c_header| {
            c_header
                .to_str()
                .map(|s| s.into())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_raw_header() -> io::Result<()> {
        const NUL: u8 = 0x00;

        let raw_header = "##fileformat=VCFv4.3\n";

        let mut data = 22u32.to_le_bytes().to_vec(); // l_text = 22
        data.extend_from_slice(raw_header.as_bytes());
        data.push(NUL);

        let mut reader = &data[..];
        let actual = read_raw_header(&mut reader)?;

        assert_eq!(actual, raw_header);

        Ok(())
    }
}
