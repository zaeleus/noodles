mod magic_number;

use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) use self::magic_number::read_magic_number;

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<vcf::Header>
where
    R: AsyncRead + Unpin,
{
    let raw_header = read_raw_header(reader).await?;

    let mut header: vcf::Header = raw_header
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *header.string_maps_mut() = raw_header
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(header)
}

async fn read_raw_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncRead + Unpin,
{
    let l_text = reader.read_u32_le().await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; l_text];
    reader.read_exact(&mut buf).await?;

    c_str_to_string(&buf)
}

fn c_str_to_string(buf: &[u8]) -> io::Result<String> {
    use std::ffi::CStr;

    CStr::from_bytes_with_nul(buf)
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

    #[tokio::test]
    async fn test_read_raw_header() -> io::Result<()> {
        let data = [
            0x08, 0x00, 0x00, 0x00, // l_text = 8
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73, 0x00, // text = b"noodles\x00"
        ];

        let mut reader = &data[..];
        assert_eq!(read_raw_header(&mut reader).await?, "noodles");

        Ok(())
    }
}
