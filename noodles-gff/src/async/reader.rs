mod lazy_line;

use futures::{stream, Stream, TryStreamExt};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use self::lazy_line::read_lazy_line;
use crate::{lazy, Directive, Line, Record};

/// An async GFF reader.
pub struct Reader<R> {
    inner: R,
    buf: String,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Creates an async GFF reader.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: String::new(),
        }
    }

    /// Reads a raw GFF line.
    pub async fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf).await
    }

    /// Reads a lazy line.
    pub async fn read_lazy_line(&mut self, line: &mut lazy::Line) -> io::Result<usize> {
        read_lazy_line(&mut self.inner, &mut self.buf, line).await
    }

    /// Returns an stream over lines.
    pub fn lines(&mut self) -> impl Stream<Item = io::Result<Line>> + '_ {
        Box::pin(stream::try_unfold(
            (self, String::new()),
            |(reader, mut buf)| async {
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

    /// Returns an stream over records.
    pub fn records(&mut self) -> impl Stream<Item = io::Result<Record>> + '_ {
        Box::pin(stream::try_unfold(self.lines(), |mut lines| async {
            loop {
                match lines.try_next().await? {
                    None | Some(Line::Directive(Directive::StartOfFasta)) => return Ok(None),
                    Some(Line::Record(record)) => return Ok(Some((record, lines))),
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
