use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_vcf as vcf;

pub(super) fn write_header<W>(writer: &mut W, header: &vcf::Header) -> io::Result<()>
where
    W: Write,
{
    let raw_header = serialize_header(header)?;
    let c_raw_header =
        CString::new(raw_header).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let text = c_raw_header.as_bytes_with_nul();
    let l_text =
        u32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32::<LittleEndian>(l_text)?;
    writer.write_all(text)?;

    Ok(())
}

pub(crate) fn serialize_header(header: &vcf::Header) -> io::Result<Vec<u8>> {
    let mut writer = vcf::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    Ok(writer.into_inner())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        use noodles_vcf::header::FileFormat;

        let mut buf = Vec::new();

        let header = vcf::Header::builder()
            .set_file_format(FileFormat::new(4, 5))
            .build();

        write_header(&mut buf, &header)?;

        let mut expected = 61i32.to_le_bytes().to_vec();

        let text = b"##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\0";
        expected.extend_from_slice(text);

        assert_eq!(buf, expected);

        Ok(())
    }
}
