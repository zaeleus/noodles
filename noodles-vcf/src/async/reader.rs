mod query;

use futures::{stream, Stream};
use memchr::memchr;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_tabix as tabix;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncSeek};

use self::query::query;
use crate::{reader::resolve_region, Header, Record};

const LINE_FEED: char = '\n';
const CARRIAGE_RETURN: char = '\r';

const HEADER_PREFIX: u8 = b'#';

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
/// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
/// use futures::TryStreamExt;
/// use noodles_vcf as vcf;
/// use tokio::{fs::File, io::BufReader};
///
/// let mut reader = File::open("sample.vcf")
///     .await
///     .map(BufReader::new)
///     .map(vcf::AsyncReader::new)?;
///
/// let header = reader.read_header().await?.parse()?;
///
/// let mut records = reader.records(&header);
///
/// while let Some(record) = records.try_next().await? {
///     // ...
/// }
/// # Ok(())
/// # }
/// ```
pub struct Reader<R> {
    inner: R,
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
    /// let reader = vcf::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the raw VCF header.
    ///
    /// This reads all header lines prefixed with a `#` (number sign), which includes the header
    /// header (`#CHROM`...).
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw VCF header as a [`String`], and as such, it is not necessarily valid.
    /// The raw header can subsequently be parsed as a [`crate::Header`].
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
    /// let mut reader = vcf::AsyncReader::new(&data[..]);
    /// let header = reader.read_header().await?;
    ///
    /// assert_eq!(header, "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<String> {
        read_header(&mut self.inner).await
    }

    /// Reads a single raw VCF record.
    ///
    /// This reads from the underlying stream until a newline is reached and appends it to the
    /// given buffer, sans the final newline. The buffer does not necessarily represent a valid VCF
    /// record but can subsequently be parsed as a [`crate::Record`].
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`] and
    /// [`Self::query`]), but using this method allows control of the line buffer and whether the
    /// raw record should be parsed.
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
    /// let mut reader = vcf::AsyncReader::new(&data[..]);
    /// reader.read_header().await?;
    ///
    /// let mut buf = String::new();
    /// reader.read_record(&mut buf).await?;
    ///
    /// assert_eq!(buf, "sq0\t1\t.\tA\t.\t.\tPASS\t.");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf).await
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
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_vcf as vcf;
    ///
    /// let data = b"##fileformat=VCFv4.3
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// sq0\t1\t.\tA\t.\t.\tPASS\t.
    /// ";
    ///
    /// let mut reader = vcf::AsyncReader::new(&data[..]);
    /// let header = reader.read_header().await?.parse()?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
    ) -> impl Stream<Item = io::Result<Record>> + 'r {
        Box::pin(stream::try_unfold(
            (&mut self.inner, String::new()),
            |(mut reader, mut buf)| async {
                buf.clear();

                match read_line(&mut reader, &mut buf).await? {
                    0 => Ok(None),
                    _ => Record::try_from_str(&buf, header)
                        .map(|record| Some((record, (reader, buf))))
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                }
            },
        ))
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead,
{
    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// use noodles_vcf as vcf;
    ///
    /// let data = [];
    /// let reader = vcf::AsyncReader::new(bgzf::AsyncReader::new(&data[..]));
    ///
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::default());
    /// # Ok::<(), io::Error>(())
    /// ```
    ///
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Seeks the underlying BGZF stream to the given virtual position.
    ///
    /// Virtual positions typically come from an associated index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bgzf as bgzf;
    /// use noodles_vcf as vcf;
    ///
    /// let data = Cursor::new([]);
    /// let mut reader = vcf::AsyncReader::new(bgzf::AsyncReader::new(data));
    ///
    /// let virtual_position = bgzf::VirtualPosition::default();
    /// reader.seek_virtual_position(virtual_position).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn seek_virtual_position(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek_virtual_position(pos).await
    }

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
    ///     .map(bgzf::AsyncReader::new)
    ///     .map(vcf::AsyncReader::new)?;
    ///
    /// let header = reader.read_header().await?.parse()?;
    ///
    /// let index = tabix::read("sample.vcf.gz.tbi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'r>(
        &'r mut self,
        header: &'r Header,
        index: &tabix::Index,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + 'r> {
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

async fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncBufRead + Unpin,
{
    let mut header_buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let buf = reader.fill_buf().await?;

        if (i == 0 || is_eol) && buf.first().map(|&b| b != HEADER_PREFIX).unwrap_or(true) {
            break;
        }

        let (read_eol, len) = if let Some(i) = memchr(LINE_FEED as u8, buf) {
            header_buf.extend(&buf[..=i]);
            (true, i + 1)
        } else {
            header_buf.extend(buf);
            (false, buf.len())
        };

        is_eol = read_eol;

        reader.consume(len);
    }

    String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
