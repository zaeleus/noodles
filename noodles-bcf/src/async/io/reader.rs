pub mod header;
mod query;
mod record;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead, AsyncSeek};

use self::{header::read_header, query::query, record::read_record};
use crate::Record;

/// An async BCF reader.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// #
/// # #[tokio::main]
/// # async fn main() -> io::Result<()> {
/// use futures::TryStreamExt;
/// use noodles_bcf as bcf;
/// use tokio::fs::File;
///
/// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
/// reader.read_header().await?;
///
/// let mut records = reader.records();
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
    R: AsyncRead + Unpin,
{
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use tokio::io;
    /// let reader = bcf::r#async::io::Reader::from(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use tokio::io;
    /// let mut reader = bcf::r#async::io::Reader::from(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use tokio::io;
    /// let reader = bcf::r#async::io::Reader::from(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Returns a BCF header reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bcf as bcf;
    /// use tokio::{fs::File, io::AsyncReadExt};
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number().await?;
    /// header_reader.read_format_version().await?;
    ///
    /// let mut raw_vcf_header_reader = header_reader.raw_vcf_header_reader().await?;
    /// let mut raw_header = String::new();
    /// raw_vcf_header_reader.read_to_string(&mut raw_header).await?;
    /// raw_vcf_header_reader.discard_to_end().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the VCF header.
    ///
    /// The BCF magic number is checked, and the file format version is discarded.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bcf as bcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<vcf::Header> {
        read_header(&mut self.inner).await
    }

    /// Reads a single record without decoding (most of) its fields.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`]), but using
    /// this method directly allows the reuse of a single [`Record`] buffer.
    ///
    /// If successful, the record size is returned. If a record size of 0 is returned, the stream
    /// reached EOF.
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bcf as bcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
    ///
    /// let mut record = bcf::Record::default();
    /// reader.read_record(&mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record).await
    }

    /// Returns an (async) stream over lazy records starting from the current (input) stream
    /// position.
    ///
    /// The (input) stream is expected to be directly after the header or at the start of another
    /// record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_bcf as bcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
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
        Box::pin(stream::try_unfold(
            (&mut self.inner, Record::default()),
            |(mut reader, mut record)| async {
                read_record(&mut reader, &mut record)
                    .await
                    .map(|n| match n {
                        0 => None,
                        _ => Some((record.clone(), (reader, record))),
                    })
            },
        ))
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async BCF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use tokio::io;
    /// let reader = bcf::r#async::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self::from(bgzf::AsyncReader::new(inner))
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Returns a stream over records that intersect the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_bcf as bcf;
    /// use noodles_core::Region;
    /// use noodles_csi as csi;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let index = csi::r#async::read("sample.bcf.csi").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<I>(
        &mut self,
        header: &vcf::Header,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + '_>
    where
        I: BinningIndex,
    {
        use crate::io::reader::resolve_region;

        let reference_sequence_id = resolve_region(header.string_maps().contigs(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(query(
            self,
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
    }
}

impl<R> From<R> for Reader<R> {
    fn from(inner: R) -> Self {
        Self { inner }
    }
}
