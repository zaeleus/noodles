mod lazy_line;

use futures::{stream, Stream, TryStreamExt};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use self::lazy_line::read_lazy_line;
use crate::{lazy, DirectiveBuf, LineBuf, Record};

/// An async GFF reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// use tokio::io;
    /// let reader = gff::r#async::io::Reader::new(io::empty());
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
    /// use noodles_gff as gff;
    /// use tokio::io;
    /// let mut reader = gff::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// use tokio::io;
    /// let reader = gff::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Creates an async GFF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// use tokio::io;
    /// let reader = gff::r#async::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a GFF line buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_gff as gff;
    ///
    /// let data = b"##gff-version 3\n";
    /// let mut reader = gff::r#async::io::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_line(&mut buf).await?;
    /// assert_eq!(buf, "##gff-version 3");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf).await
    }

    /// Reads a lazy line.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_gff::{self as gff, lazy};
    ///
    /// let data = b"##gff-version 3\n";
    /// let mut reader = gff::r#async::io::Reader::new(&data[..]);
    ///
    /// let mut line = lazy::Line::default();
    ///
    /// reader.read_lazy_line(&mut line).await?;
    /// assert_eq!(line.kind(), lazy::line::Kind::Directive);
    ///
    /// assert_eq!(reader.read_lazy_line(&mut line).await?, 0);
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_lazy_line(&mut self, line: &mut lazy::Line) -> io::Result<usize> {
        read_lazy_line(&mut self.inner, line).await
    }

    /// Returns a stream over line buffers.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_gff::{self as gff, LineBuf};
    ///
    /// let data = b"##gff-version 3\n";
    /// let mut reader = gff::r#async::io::Reader::new(&data[..]);
    /// let mut lines = reader.lines();
    ///
    /// let line = lines.try_next().await?;
    /// assert!(matches!(line, Some(LineBuf::Directive(_))));
    ///
    /// assert!(lines.try_next().await?.is_none());
    /// # Ok(())
    /// # }
    /// ```
    pub fn lines(&mut self) -> impl Stream<Item = io::Result<LineBuf>> + '_ {
        Box::pin(stream::try_unfold(
            (self, String::new()),
            |(reader, mut buf)| async {
                buf.clear();

                reader.read_line(&mut buf).await.and_then(|n| match n {
                    0 => Ok(None),
                    _ => match buf.parse() {
                        Ok(line) => Ok(Some((line, (reader, buf)))),
                        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
                    },
                })
            },
        ))
    }

    /// Returns a stream over records.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_gff as gff;
    /// use futures::TryStreamExt;
    ///
    /// let data = b"##gff-version 3\n";
    /// let mut reader = gff::r#async::io::Reader::new(&data[..]);
    /// let mut records = reader.records();
    ///
    /// assert!(records.try_next().await?.is_none());
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(&mut self) -> impl Stream<Item = io::Result<Record>> + '_ {
        Box::pin(stream::try_unfold(self.lines(), |mut lines| async {
            loop {
                match lines.try_next().await? {
                    None | Some(LineBuf::Directive(DirectiveBuf::StartOfFasta)) => return Ok(None),
                    Some(LineBuf::Record(record)) => return Ok(Some((record, lines))),
                    _ => {}
                }
            }
        }))
    }
}

async fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

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
