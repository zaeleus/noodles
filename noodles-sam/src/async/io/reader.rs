mod header;

use futures::{stream, Stream};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use self::header::read_header;
use crate::{alignment::RecordBuf, Header};

/// An async SAM reader.
pub struct Reader<R> {
    inner: R,
    buf: Vec<u8>,
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

    /// Reads a single SAM record.
    ///
    /// This reads a line from the underlying stream until a newline is reached and parses that
    /// line into the given record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`]), but using this
    /// method allows control of the line buffer and whether the raw record should be parsed.
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
    /// reader.read_record(&header, &mut record).await?;
    ///
    /// assert_eq!(record, RecordBuf::default());
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(
        &mut self,
        header: &Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        read_record(&mut self.inner, &mut self.buf, header, record).await
    }

    /// Returns an (async) stream over records starting from the current (input) stream position.
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
    /// let mut records = reader.records(&header);
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> impl Stream<Item = io::Result<RecordBuf>> + 'a {
        Box::pin(stream::try_unfold(
            (self, RecordBuf::default()),
            |(reader, mut record)| async {
                match reader.read_record(header, &mut record).await? {
                    0 => Ok(None),
                    _ => Ok(Some((record.clone(), (reader, record)))),
                }
            },
        ))
    }
}

async fn read_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    header: &Header,
    record: &mut RecordBuf,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    use crate::io::reader::record_buf::parse_record_buf;

    buf.clear();

    match read_line(reader, buf).await? {
        0 => Ok(0),
        n => {
            parse_record_buf(buf, header, record)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
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
