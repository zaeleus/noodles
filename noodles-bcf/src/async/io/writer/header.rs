use std::ffi::CString;

use noodles_vcf as vcf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub(super) async fn write_header<W>(writer: &mut W, header: &vcf::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::io::writer::header::serialize_header;

    let raw_header = serialize_header(header)?;
    let c_raw_header =
        CString::new(raw_header).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let text = c_raw_header.as_bytes_with_nul();
    let l_text =
        u32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    writer.write_u32_le(l_text).await?;
    writer.write_all(text).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();
        let header = vcf::Header::default();
        write_header(&mut buf, &header).await?;

        let mut expected = 61i32.to_le_bytes().to_vec();

        let text = b"##fileformat=VCFv4.4\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\0";
        expected.extend_from_slice(text);

        assert_eq!(buf, expected);

        Ok(())
    }
}
