pub mod header;
mod query;
mod record;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncSeek};

use self::{header::read_header, query::query, record::read_record};
use crate::{io::reader::resolve_region, variant::RecordBuf, Header, Record};

const LINE_FEED: char = '\n';
const CARRIAGE_RETURN: char = '\r';

/// An async VCF reader.
///
/// The VCF format has two main parts: 1) a header and 2) a list of VCF records.
///
/// Each header line is prefixed with a `#` (number sign) and is terminated by the header header
/// (`#CHROM`...; inclusive).
///
/// VCF records are line-based and follow directly after the header until EOF.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> std::io::Result<()> {
/// use futures::TryStreamExt;
/// use noodles_vcf as vcf;
/// use tokio::{fs::File, io::BufReader};
///
/// let mut reader = File::open("sample.vcf")
///     .await
///     .map(BufReader::new)
///     .map(vcf::r#async::io::Reader::new)?;
///
/// let _header = reader.read_header().await?;
///
/// let mut records = reader.records();
///
/// while let Some(_record) = records.try_next().await? {
///     // ...
/// }
/// # Ok(())
/// # }
/// ```
pub struct Reader<R> {
    inner: R,
    buf: String,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let reader = vcf::r#async::io::Reader::new(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let mut reader = vcf::r#async::io::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let reader = vcf::r#async::io::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Creates an async VCF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let data = [];
    /// let reader = vcf::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: String::new(),
        }
    }

    /// Returns an async VCF header reader.
    ///
    /// This creates an adapter that reads at most the length of the header, i.e., all lines
    /// prefixed with a `#` (number sign).
    ///
    /// It is more ergonomic to read and parse the header using [`Self::read_header`], but using
    /// this adapter allows for control of how the header is read, e.g., to read the raw VCF
    /// header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_vcf as vcf;
    /// use tokio::io::AsyncReadExt;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::r#async::io::Reader::new(&data[..]);
    /// let mut header_reader = reader.header_reader();
    ///
    /// let mut raw_header = String::new();
    /// header_reader.read_to_string(&mut raw_header).await?;
    ///
    /// assert_eq!(raw_header, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    /// # Ok(())
    /// # }
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the VCF header.
    ///
    /// This reads all header lines prefixed with a `#` (number sign), which includes the header
    /// header (`#CHROM`...), and parses it as a [`crate::Header`].
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::r#async::io::Reader::new(&data[..]);
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<Header> {
        read_header(&mut self.inner).await
    }

    /// Reads a single VCF record.
    ///
    /// This reads a line from the underlying stream until a newline is reached and parses that
    /// line into the given record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`] and
    /// [`Self::query`]), but using this method allows control of the record buffer.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::r#async::io::Reader::new(&data[..]);
    /// let header = reader.read_header().await?;
    ///
    /// let mut record = vcf::variant::RecordBuf::default();
    /// reader.read_record_buf(&header, &mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record_buf(
        &mut self,
        header: &Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        use crate::io::reader::parse_record_buf;

        self.buf.clear();

        match read_line(&mut self.inner, &mut self.buf).await? {
            0 => Ok(0),
            n => {
                parse_record_buf(&self.buf, header, record)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                Ok(n)
            }
        }
    }

    /// Reads a single record without eagerly parsing its fields.
    ///
    /// The reads VCF record fields from the underlying stream into the given record's buffer until
    /// a newline is reached. No fields are parsed, meaning the record is not necessarily valid.
    /// However, the structure of the line is guaranteed to be record-like.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::r#async::io::Reader::new(&data[..]);
    /// reader.read_header().await?;
    ///
    /// let mut record = vcf::Record::default();
    /// reader.read_record(&mut record).await?;
    /// # Ok::<_, std::io::Error>(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, &mut self.buf, record).await
    }

    /// Returns a stream over records.
    ///
    /// The (input) stream is expected to be directly after the header or at the start of another
    /// record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_vcf as vcf;
    ///
    /// const DATA: &[u8] = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::r#async::io::Reader::new(DATA);
    /// let header = reader.read_header().await?;
    ///
    /// let mut records = reader.records();
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(&mut self) -> impl Stream<Item = io::Result<Record>> + '_ {
        Box::pin(stream::try_unfold(self, move |reader| async move {
            let mut record = Record::default();

            reader.read_record(&mut record).await.map(|n| match n {
                0 => None,
                _ => Some((record, reader)),
            })
        }))
    }

    /// Returns an (async) stream over records starting from the current (input) stream position.
    ///
    /// The (input) stream is expected to be directly after the header or at the start of another
    /// record.
    ///
    /// Unlike [`Self::read_record`], each record is parsed as a [`Record`].
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::r#async::io::Reader::new(&data[..]);
    /// let header = reader.read_header().await?;
    ///
    /// let mut records = reader.record_bufs(&header);
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn record_bufs<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<RecordBuf>> + 'r {
        Box::pin(stream::try_unfold(self, move |reader| async move {
            let mut record = RecordBuf::default();

            reader
                .read_record_buf(header, &mut record)
                .await
                .map(|n| match n {
                    0 => None,
                    _ => Some((record, reader)),
                })
        }))
    }
}

impl<R> Reader<bgzf::r#async::io::Reader<R>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Returns a stream over records that intersects the given region.
    ///
    /// The position of the (input) stream is expected to after the header or at the start of
    /// another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_bgzf as bgzf;
    /// use noodles_core::Region;
    /// use noodles_tabix as tabix;
    /// use noodles_vcf as vcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.vcf.gz")
    ///     .await
    ///     .map(bgzf::r#async::io::Reader::new)
    ///     .map(vcf::r#async::io::Reader::new)?;
    ///
    /// let header = reader.read_header().await?;
    ///
    /// let index = tabix::r#async::fs::read("sample.vcf.gz.tbi").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'r, I>(
        &'r mut self,
        header: &'r Header,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + 'r>
    where
        I: BinningIndex,
    {
        let (reference_sequence_id, reference_sequence_name) = resolve_region(index, region)?;

        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(query(
            self,
            chunks,
            reference_sequence_name,
            region.interval(),
            header,
        ))
    }
}

async fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    match reader.read_line(buf).await {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
                    buf.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_line() -> io::Result<()> {
        async fn t(buf: &mut String, mut data: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut data, buf).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = String::new();

        t(&mut buf, b"noodles\n", "noodles").await?;
        t(&mut buf, b"noodles\r\n", "noodles").await?;
        t(&mut buf, b"noodles", "noodles").await?;

        Ok(())
    }
}
