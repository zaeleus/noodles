mod format_version;
mod magic_number;
mod vcf_header;

use noodles_vcf::{self as vcf, header::StringMaps};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt};

use self::{format_version::read_format_version, magic_number::read_magic_number};
use crate::MAGIC_NUMBER;

struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    fn new(inner: R) -> Self {
        Self { inner }
    }

    async fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        read_magic_number(&mut self.inner).await
    }

    async fn read_format_version(&mut self) -> io::Result<(u8, u8)> {
        read_format_version(&mut self.inner).await
    }

    async fn raw_vcf_header_reader(&mut self) -> io::Result<vcf_header::Reader<&mut R>> {
        let len = self.inner.read_u32_le().await.map(u64::from)?;
        Ok(vcf_header::Reader::new(&mut self.inner, len))
    }
}

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<vcf::Header>
where
    R: AsyncRead + Unpin,
{
    let mut header_reader = Reader::new(reader);
    read_header_inner(&mut header_reader).await
}

async fn read_header_inner<R>(reader: &mut Reader<R>) -> io::Result<vcf::Header>
where
    R: AsyncRead + Unpin,
{
    reader
        .read_magic_number()
        .await
        .and_then(crate::io::reader::header::magic_number::validate)?;

    reader.read_format_version().await?;

    let mut raw_vcf_header_reader = reader.raw_vcf_header_reader().await?;
    read_vcf_header(&mut raw_vcf_header_reader).await
}

async fn read_vcf_header<R>(reader: &mut vcf_header::Reader<R>) -> io::Result<vcf::Header>
where
    R: AsyncRead + Unpin,
{
    let mut parser = vcf::header::Parser::default();
    let mut string_maps = StringMaps::default();

    let mut buf = Vec::new();

    while read_line(reader, &mut buf).await? != 0 {
        let entry = parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        string_maps
            .insert_entry(&entry)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    reader.discard_to_end().await?;

    let mut header = parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *header.string_maps_mut() = string_maps;

    Ok(header)
}

async fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst).await? {
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

    #[tokio::test]
    async fn test_read_header() -> io::Result<()> {
        use vcf::header::FileFormat;

        const NUL: u8 = 0x00;

        let raw_header = b"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
";

        let mut data = vec![
            b'B', b'C', b'F', // magic
            0x02, 0x02, // major_version, minor_version
        ];
        data.extend(61u32.to_le_bytes()); // l_text
        data.extend(raw_header); // text
        data.push(NUL);

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        let expected = vcf::Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
