use futures::{Stream, stream};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncReadExt};

use crate::Record;

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

/// An async FASTQ reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let reader = fastq::r#async::io::Reader::new(io::empty());
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
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let mut reader = fastq::r#async::io::Reader::new(io::empty());
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
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let reader = fastq::r#async::io::Reader::new(io::empty());
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
    /// Creates an async FASTQ reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// use tokio::io;
    /// let reader = fastq::r#async::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a FASTQ record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_fastq as fastq;
    ///
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let mut reader = fastq::r#async::io::Reader::new(&data[..]);
    ///
    /// let mut record = fastq::Record::default();
    /// reader.read_record(&mut record).await?;
    ///
    /// assert_eq!(record.name(), &b"r0"[..]);
    /// assert_eq!(record.sequence(), b"ATCG");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record).await
    }

    /// Returns an (async) stream over records starting from the current (input) stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_fastq as fastq;
    ///
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let mut reader = fastq::r#async::io::Reader::new(&data[..]);
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
            |(mut reader, mut buf)| async {
                read_record(&mut reader, &mut buf).await.map(|n| match n {
                    0 => None,
                    _ => Some((buf.clone(), (reader, buf))),
                })
            },
        ))
    }
}

async fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    record.clear();

    let mut len = match read_name(reader, record).await? {
        0 => return Ok(0),
        n => n,
    };

    len += read_line(reader, record.sequence_mut()).await?;
    len += read_description(reader, &mut Vec::new()).await?;
    len += read_line(reader, record.quality_scores_mut()).await?;

    Ok(len)
}

async fn read_name<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    use memchr::memchr2;

    const NAME_PREFIX: u8 = b'@';

    const HORIZONTAL_TAB: u8 = b'\t';
    const SPACE: u8 = b' ';

    match reader.read_u8().await {
        Ok(NAME_PREFIX) => {
            let n = read_line(reader, record.name_mut()).await.map(|n| n + 1)?;

            if let Some(i) = memchr2(SPACE, HORIZONTAL_TAB, record.name()) {
                let description = record.name_mut().split_off(i + 1);
                record.name_mut().pop();
                *record.description_mut() = description.into();
            }

            Ok(n)
        }
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid name prefix",
        )),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

async fn read_description<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const DESCRIPTION_PREFIX: u8 = b'+';

    match reader.read_u8().await? {
        DESCRIPTION_PREFIX => read_line(reader, buf).await.map(|n| n + 1),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid description prefix",
        )),
    }
}

async fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
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
    use crate::record::Definition;

    #[tokio::test]
    async fn test_read_record() -> io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles:2/1
TCGA
+noodles:2/1
dcba
";

        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record).await?;
        let expected = Record::new(Definition::new("noodles:1/1", ""), "AGCT", "abcd");
        assert_eq!(record, expected);

        read_record(&mut reader, &mut record).await?;
        let expected = Record::new(Definition::new("noodles:2/1", ""), "TCGA", "dcba");
        assert_eq!(record, expected);

        let n = read_record(&mut reader, &mut record).await?;
        assert_eq!(n, 0);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_name() -> io::Result<()> {
        let mut record = Record::default();

        let data = b"@r0\n";
        let mut reader = &data[..];
        record.clear();
        read_name(&mut reader, &mut record).await?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert!(record.description().is_empty());

        let data = b"@r0 LN:4\n";
        let mut reader = &data[..];
        record.clear();
        read_name(&mut reader, &mut record).await?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.description(), &b"LN:4"[..]);

        let data = b"@r0\tLN:4\n";
        let mut reader = &data[..];
        record.clear();
        read_name(&mut reader, &mut record).await?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.description(), &b"LN:4"[..]);

        let data = b"r0\n";
        let mut reader = &data[..];
        record.clear();
        assert!(matches!(
            read_name(&mut reader, &mut record).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_description() -> io::Result<()> {
        let mut buf = Vec::new();

        let data = b"+r0\n";
        let mut reader = &data[..];
        buf.clear();
        read_description(&mut reader, &mut buf).await?;
        assert_eq!(buf, b"r0");

        let data = b"r0\n";
        let mut reader = &data[..];
        buf.clear();
        assert!(matches!(
            read_description(&mut reader, &mut buf).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

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
