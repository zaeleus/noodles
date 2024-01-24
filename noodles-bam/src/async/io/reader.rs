mod header;
mod query;
mod record_buf;

use bytes::BytesMut;
use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::{self as sam, alignment::RecordBuf};
use tokio::io::{self, AsyncRead, AsyncReadExt, AsyncSeek};

use self::{
    header::read_header,
    query::query,
    record_buf::{read_block_size, read_record_buf},
};
use crate::{io::reader::resolve_region, Record, MAGIC_NUMBER};

/// An async BAM reader.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> std::io::Result<()> {
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
pub struct Reader<R> {
    inner: R,
    buf: BytesMut,
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::r#async::io::Reader::from(&data[..]);
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let mut reader = bam::r#async::io::Reader::from(&data[..]);
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::r#async::io::Reader::from(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
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
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<sam::Header> {
        read_magic(&mut self.inner).await?;
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
    /// # async fn main() -> std::io::Result<()> {
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
        header: &sam::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, header, &mut self.buf, record).await
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
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
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
        let block_size = match read_block_size(&mut self.inner).await? {
            0 => return Ok(0),
            n => n,
        };

        let fields = record.fields_mut();

        fields.buf.resize(block_size, 0);
        self.inner.read_exact(&mut fields.buf).await?;

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
    /// # async fn main() -> std::io::Result<()> {
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
    pub fn record_bufs<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> impl Stream<Item = io::Result<RecordBuf>> + '_ {
        Box::pin(stream::try_unfold(
            (&mut self.inner, &mut self.buf, RecordBuf::default()),
            move |(reader, buf, mut record)| async move {
                read_record_buf(reader, header, buf, &mut record)
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
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
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
    pub fn records(&mut self) -> impl Stream<Item = io::Result<Record>> + '_ {
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

impl<R> Reader<bgzf::AsyncReader<R>>
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
        Self::from(bgzf::AsyncReader::new(reader))
    }

    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_bgzf as bgzf;
    ///
    /// let data = Vec::new();
    /// let reader = bam::r#async::io::Reader::new(&data[..]);
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Seeks the underlying BGZF reader to the given virtual position.
    ///
    /// Virtual positions typically come from the associated BAM index file.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_bgzf as bgzf;
    ///
    /// let data = [];
    /// let mut reader = bam::r#async::io::Reader::new(Cursor::new(data));
    ///
    /// let virtual_position = bgzf::VirtualPosition::default();
    /// reader.seek(virtual_position).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos).await
    }

    /// Returns a stream over records that intersect the given region.
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
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let index = bai::r#async::read("sample.bam.bai").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'a, I>(
        &'a mut self,
        header: &'a sam::Header,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<RecordBuf>> + '_>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(header.reference_sequences(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(query(
            self,
            header,
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
    }
}

impl<R> From<R> for Reader<R> {
    fn from(inner: R) -> Self {
        Self {
            inner,
            buf: BytesMut::new(),
        }
    }
}

async fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic() {
        let data = b"BAM\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).await.is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
