mod header;
mod query;
mod record;
mod record_buf;

use futures::{stream, Stream, TryStreamExt};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncSeek};

use self::{header::read_header, query::query, record::read_record, record_buf::read_record_buf};
use crate::{alignment::RecordBuf, Header, Record};

/// An async SAM reader.
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
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::r#async::io::Reader::new(&data[..]);
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
    /// use noodles_sam as sam;
    /// let data = [];
    /// let mut reader = sam::r#async::io::Reader::new(&data[..]);
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
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::r#async::io::Reader::new(&data[..]);
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
    /// Creates an async SAM reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }

    /// Reads the SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// The SAM header is optional, and if it is missing, an empty [`Header`] is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::r#async::io::Reader::new(&data[..]);
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<Header> {
        read_header(&mut self.inner).await
    }

    /// Reads a record into an alignment record buffer.
    ///
    /// This reads a line from the underlying stream until a newline is reached and parses that
    /// line into the given record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::records`]), but using
    /// this method directly allows reuse of a [`RecordBuf`].
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::r#async::io::Reader::new(&data[..]);
    /// let header = reader.read_header().await?;
    ///
    /// let mut record = RecordBuf::default();
    /// reader.read_record_buf(&header, &mut record).await?;
    ///
    /// assert_eq!(record, RecordBuf::default());
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record_buf(
        &mut self,
        header: &Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record_buf(&mut self.inner, &mut self.buf, header, record).await
    }

    /// Returns an (async) stream over alignment record buffers starting from the current (input)
    /// stream position.
    ///
    /// The (input) stream is expected to be directly after the header or at the start of another
    /// record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::r#async::io::Reader::new(&data[..]);
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
        header: &'a Header,
    ) -> impl Stream<Item = io::Result<RecordBuf>> + 'a {
        Box::pin(stream::try_unfold(
            (self, RecordBuf::default()),
            |(reader, mut record)| async {
                match reader.read_record_buf(header, &mut record).await? {
                    0 => Ok(None),
                    _ => Ok(Some((record.clone(), (reader, record)))),
                }
            },
        ))
    }

    /// Reads a record.
    ///
    /// This reads SAM fields from the underlying stream into the given record's buffer until a
    /// newline is reached. No fields are decoded, meaning the record is not necessarily valid.
    /// However, the structure of the buffer is guaranteed to be record-like.
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
    /// use noodles_sam as sam;
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::r#async::io::Reader::new(&data[..]);
    /// let header = reader.read_header().await?;
    ///
    /// let mut record = sam::Record::default();
    /// reader.read_record(&mut record).await?;
    ///
    /// assert_eq!(record, sam::Record::default());
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, &mut self.buf, record).await
    }

    /// Returns an (async) stream over records.
    ///
    /// The (input) stream is expected to be directly after the header or at the start of a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use tokio::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = sam::r#async::io::Reader::new(io::empty());
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
        Box::pin(stream::try_unfold(
            (self, Record::default()),
            |(reader, mut record)| async {
                match reader.read_record(&mut record).await? {
                    0 => Ok(None),
                    _ => Ok(Some((record.clone(), (reader, record)))),
                }
            },
        ))
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
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
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi as csi;
    /// use noodles_sam as sam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.sam")
    ///     .await
    ///     .map(bgzf::AsyncReader::new)
    ///     .map(sam::r#async::io::Reader::new)?;
    ///
    /// let header = reader.read_header().await?;
    ///
    /// let index = csi::r#async::read("sample.sam.csi").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'r, 'h: 'r, I>(
        &'r mut self,
        header: &'h Header,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + 'r>
    where
        I: BinningIndex,
    {
        use crate::io::reader::resolve_region;

        let reference_sequence_id = resolve_region(header.reference_sequences(), region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(query(
            self,
            chunks,
            header,
            reference_sequence_id,
            region.interval(),
        ))
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    ///
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi as csi;
    /// use noodles_sam as sam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.sam")
    ///     .await
    ///     .map(bgzf::AsyncReader::new)
    ///     .map(sam::r#async::io::Reader::new)?;
    ///
    /// let index = csi::r#async::read("sample.sam.csi").await?;
    /// let mut query = reader.query_unmapped(&index).await?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn query_unmapped<I>(
        &mut self,
        index: &I,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + '_>
    where
        I: BinningIndex,
    {
        if let Some(pos) = index.last_first_record_start_position() {
            self.get_mut().seek(pos).await?;
        } else {
            self.seek_to_first_record().await?;
        }

        Ok(Box::pin(stream::try_unfold(
            self.records(),
            |mut records| async {
                loop {
                    match records.try_next().await? {
                        Some(record) => {
                            if record.flags()?.is_unmapped() {
                                return Ok(Some((record, records)));
                            }
                        }
                        None => return Ok(None),
                    }
                }
            },
        )))
    }
}

async fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    match reader.read_until(LINE_FEED, buf).await? {
        0 => Ok(0),
        n => {
            if buf.ends_with(&[LINE_FEED]) {
                buf.pop();

                if buf.ends_with(&[CARRIAGE_RETURN]) {
                    buf.pop();
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
    async fn test_read_line() -> io::Result<()> {
        async fn t(buf: &mut Vec<u8>, mut data: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_line(&mut data, buf).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"noodles\n", b"noodles").await?;
        t(&mut buf, b"noodles\r\n", b"noodles").await?;
        t(&mut buf, b"noodles", b"noodles").await?;

        Ok(())
    }
}
