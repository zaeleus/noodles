use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::Record;

/// An async FASTQ writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let writer = fastq::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let mut writer = fastq::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let writer = fastq::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async FASTQ writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let writer = fastq::r#async::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a FASTQ record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_fastq::{self as fastq, record::Definition};
    ///
    /// let mut writer = fastq::r#async::io::Writer::new(Vec::new());
    ///
    /// let record = fastq::Record::new(Definition::new("r0", ""), "ATCG", "NDLS");
    /// writer.write_record(&record).await?;
    ///
    /// assert_eq!(writer.get_ref(), b"@r0\nATCG\n+\nNDLS\n");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, record).await
    }
}

async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const NAME_PREFIX: &[u8] = b"@";
    const DEFINITION_SEPARATOR: &[u8] = b" ";
    const LINE_FEED: &[u8] = b"\n";

    writer.write_all(NAME_PREFIX).await?;
    writer.write_all(record.name()).await?;

    if !record.description().is_empty() {
        writer.write_all(DEFINITION_SEPARATOR).await?;
        writer.write_all(record.description()).await?;
    }

    writer.write_all(LINE_FEED).await?;

    writer.write_all(record.sequence()).await?;
    writer.write_all(LINE_FEED).await?;

    writer.write_all(b"+").await?;
    writer.write_all(LINE_FEED).await?;

    writer.write_all(record.quality_scores()).await?;
    writer.write_all(LINE_FEED).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Definition;

    #[tokio::test]
    async fn test_write_record() -> io::Result<()> {
        let mut buf = Vec::new();

        let mut record = Record::new(Definition::new("r0", ""), "ACGT", "NDLS");
        write_record(&mut buf, &record).await?;
        let expected = b"@r0\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        record.description_mut().extend_from_slice(b"LN:4");

        buf.clear();
        write_record(&mut buf, &record).await?;
        let expected = b"@r0 LN:4\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
