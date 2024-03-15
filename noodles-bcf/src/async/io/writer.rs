mod header;

use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, header::StringMaps};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::header::write_header;
use crate::Record;

/// An async BCF writer.
pub struct Writer<W> {
    inner: W,
    string_maps: StringMaps,
    buf: Vec<u8>,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Writes a VCF header.
    pub async fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        write_file_format(&mut self.inner).await?;

        self.string_maps = StringMaps::try_from(header)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        write_header(&mut self.inner, header).await
    }

    /// Writes a record.
    pub async fn write_record(&mut self, header: &vcf::Header, record: &Record) -> io::Result<()> {
        self.write_variant_record(header, record).await
    }

    /// Writes a variant record.
    pub async fn write_variant_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> io::Result<()> {
        use crate::io::writer::write_record;

        self.buf.clear();
        write_record(&mut self.buf, header, &self.string_maps, record)?;
        self.inner.write_all(&self.buf).await
    }
}

impl<W> Writer<bgzf::AsyncWriter<W>>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BCF writer.
    ///
    /// The given stream is wrapped in a BGZF encoder.
    pub fn new(inner: W) -> Self {
        Self::from(bgzf::AsyncWriter::new(inner))
    }
}

impl<W> From<W> for Writer<W> {
    fn from(inner: W) -> Self {
        Self {
            inner,
            string_maps: StringMaps::default(),
            buf: Vec::new(),
        }
    }
}

async fn write_file_format<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::{
        io::writer::{MAJOR, MINOR},
        MAGIC_NUMBER,
    };

    writer.write_all(MAGIC_NUMBER).await?;
    writer.write_u8(MAJOR).await?;
    writer.write_u8(MINOR).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_file_format() -> io::Result<()> {
        let mut buf = Vec::new();
        write_file_format(&mut buf).await?;

        let expected = [
            b'B', b'C', b'F', // magic
            0x02, // major
            0x02, // minor
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
