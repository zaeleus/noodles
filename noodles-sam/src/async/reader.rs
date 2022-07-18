use futures::{stream, Stream};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use crate::{alignment::Record, Header};

/// An async SAM reader.
pub struct Reader<R> {
    inner: R,
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
    /// let reader = sam::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let data = [];
    /// let reader = sam::AsyncReader::new(&data[..]);
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
    /// let mut reader = sam::AsyncReader::new(&data[..]);
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
    /// let reader = sam::AsyncReader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the raw SAM header.
    ///
    /// This reads all header lines prefixed with a `@` (at sign).
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw SAM header as a [`String`], and as such, it is no necessarily valid.
    /// The raw header can subsequently be parsed as a [`crate::Header`].
    ///
    /// The SAM header is optional, and if it is missing, an empty string is returned.
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
    /// let mut reader = sam::AsyncReader::new(&data[..]);
    /// let header = reader.read_header().await?;
    ///
    /// assert_eq!(header, "@HD\tVN:1.6\n");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<String> {
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
    /// use noodles_sam::{self as sam, alignment::Record};
    ///
    /// let data = b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = sam::AsyncReader::new(&data[..]);
    /// let header = reader.read_header().await?.parse()?;
    ///
    /// let mut record = Record::default();
    /// reader.read_record(&header, &mut record).await?;
    ///
    /// assert_eq!(record, Record::default());
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, header: &Header, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, header, record).await
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
    /// let mut reader = sam::AsyncReader::new(&data[..]);
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
    pub fn records<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> impl Stream<Item = io::Result<Record>> + 'a {
        Box::pin(stream::try_unfold(
            (&mut self.inner, Record::default()),
            |(mut reader, mut record)| async {
                match read_record(&mut reader, header, &mut record).await? {
                    0 => Ok(None),
                    _ => Ok(Some((record.clone(), (reader, record)))),
                }
            },
        ))
    }
}

async fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncBufRead + Unpin,
{
    const PREFIX: u8 = b'@';
    const LINE_FEED: u8 = b'\n';

    let mut header_buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let buf = reader.fill_buf().await?;

        if (i == 0 || is_eol) && buf.first().map(|&b| b != PREFIX).unwrap_or(true) {
            break;
        }

        let (read_eol, len) = if let Some(i) = buf.iter().position(|&b| b == LINE_FEED) {
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

async fn read_record<R>(reader: &mut R, header: &Header, record: &mut Record) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    use crate::reader::record::parse_record;

    let mut buf = Vec::new();

    match read_line(reader, &mut buf).await? {
        0 => Ok(0),
        n => {
            parse_record(&buf, header, record)?;
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
