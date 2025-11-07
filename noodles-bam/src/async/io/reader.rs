pub mod header;
mod query;
mod record;
mod record_buf;

use std::future;

use futures::{Stream, StreamExt, stream};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::{self as sam, alignment::RecordBuf};
use tokio::io::{self, AsyncRead, AsyncSeek};

pub use self::query::Query;
use self::{header::read_header, record::read_record, record_buf::read_record_buf};
use crate::{Record, io::reader::resolve_region};

/// An async BAM reader.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use futures::TryStreamExt;
/// use noodles_bam as bam;
/// use tokio::fs::File;
///
/// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
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
    buf: Vec<u8>,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let reader = bam::r#async::io::Reader::from(io::empty());
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
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let mut reader = bam::r#async::io::Reader::from(io::empty());
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
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let reader = bam::r#async::io::Reader::from(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Returns a BAM header reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use tokio::{fs::File, io::AsyncReadExt};
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number().await?;
    ///
    /// let mut raw_sam_header_reader = header_reader.raw_sam_header_reader().await?;
    /// let mut raw_header = String::new();
    /// raw_sam_header_reader.read_to_string(&mut raw_header).await?;
    /// raw_sam_header_reader.discard_to_end().await?;
    ///
    /// header_reader.read_reference_sequences().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the SAM header.
    ///
    /// This verifies the BAM magic number, reads and parses the raw SAM header, and reads the
    /// binary reference sequences. If the SAM header has a reference sequence dictionary, it must
    /// match the binary reference sequences; otherwise, the binary reference sequences are added
    /// to the SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<sam::Header> {
        read_header(&mut self.inner).await
    }

    /// Reads a record into an alignment record buffer.
    ///
    /// The record block size (`bs`) is read from the underlying stream and `bs` bytes are read
    /// into an internal buffer. This buffer is then used to decode fields into the given record.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`] and
    /// [`Self::query`]), but using this method directly allows the reuse of a [`RecordBuf`].
    ///
    /// If successful, the record block size is returned. If a block size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam::alignment::RecordBuf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let mut record = RecordBuf::default();
    /// reader.read_record_buf(&header, &mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record_buf(
        &mut self,
        _header: &sam::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, &mut self.buf, record).await
    }

    /// Reads a record.
    ///
    /// The record block size (`bs`) is read from the underlying stream and `bs` bytes are read
    /// into the record's buffer. No fields are decoded, meaning the record is not necessarily
    /// valid. However, the structure of the buffer is guaranteed to be record-like.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// If successful, the record block size is returned. If a block size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
    ///
    /// let mut record = bam::Record::default();
    /// reader.read_record(&mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let fields = record.fields_mut();

        let block_size = match read_record(&mut self.inner, &mut fields.buf).await? {
            0 => return Ok(0),
            n => n,
        };

        fields.index()?;

        Ok(block_size)
    }

    /// Returns an (async) stream over alignment record buffers starting from the current (input)
    /// stream position.
    ///
    /// The (input) stream is expected to be directly after the reference sequences or at the start
    /// of another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
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
        _header: &'h sam::Header,
    ) -> impl Stream<Item = io::Result<RecordBuf>> + use<'r, R> {
        Box::pin(stream::try_unfold(
            (&mut self.inner, &mut self.buf, RecordBuf::default()),
            move |(reader, buf, mut record)| async move {
                read_record_buf(reader, buf, &mut record)
                    .await
                    .map(|n| match n {
                        0 => None,
                        _ => Some((record.clone(), (reader, buf, record))),
                    })
            },
        ))
    }

    /// Returns a stream over records.
    ///
    /// The (input) stream is expected to be directly after the reference sequences or at the start
    /// of another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
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
    pub fn records(&mut self) -> impl Stream<Item = io::Result<Record>> {
        Box::pin(stream::try_unfold(
            (self, Record::default()),
            |(this, mut record)| async {
                this.read_record(&mut record).await.map(|n| match n {
                    0 => None,
                    _ => Some((record.clone(), (this, record))),
                })
            },
        ))
    }
}

impl<R> Reader<bgzf::r#async::io::Reader<R>>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async BAM reader.
    ///
    /// The given reader must be a raw BGZF stream, as the underlying reader wraps it in a decoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(reader: R) -> Self {
        Self::from(bgzf::r#async::io::Reader::new(reader))
    }
}

impl<R> Reader<bgzf::r#async::io::Reader<R>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    // Seeks to the first record by setting the cursor to the beginning of the stream and
    // (re)reading the header.
    async fn seek_to_first_record(&mut self) -> io::Result<bgzf::VirtualPosition> {
        self.get_mut()
            .seek(bgzf::VirtualPosition::default())
            .await?;

        self.read_header().await?;

        Ok(self.get_ref().virtual_position())
    }

    /// Returns a stream over records that intersect the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_bam::{self as bam, bai};
    /// use noodles_core::Region;
    /// use noodles_sam as sam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let index = bai::r#async::fs::read("sample.bam.bai").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?.records();
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'r, I>(
        &'r mut self,
        header: &sam::Header,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'r, R>>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(header.reference_sequences(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            self.get_mut(),
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
    }

    /// Returns a stream of unmapped records after querying for the unmapped region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_bam::{self as bam, bai};
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    ///
    /// let index = bai::r#async::fs::read("sample.bam.bai").await?;
    /// let mut query = reader.query_unmapped(&index).await?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn query_unmapped<'r, I>(
        &'r mut self,
        index: &I,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + use<'r, I, R>>
    where
        I: BinningIndex,
    {
        if let Some(pos) = index.last_first_record_start_position() {
            self.get_mut().seek(pos).await?;
        } else {
            self.seek_to_first_record().await?;
        }

        Ok(self.records().filter(|result| {
            future::ready(
                result
                    .as_ref()
                    .map(|record| record.flags().is_unmapped())
                    .unwrap_or(true),
            )
        }))
    }
}

impl<R> From<R> for Reader<R> {
    fn from(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }
}
