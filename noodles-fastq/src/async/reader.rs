use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncReadExt};

use crate::Record;

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

/// An async FASTQ reader.
pub struct Reader<R> {
    inner: R,
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
    /// let data = [];
    /// let reader = fastq::AsyncReader::new(&data[..]);
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
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_fastq as fastq;
    ///
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let mut reader = fastq::AsyncReader::new(&data[..]);
    ///
    /// let mut record = fastq::Record::default();
    /// reader.read_record(&mut record).await?;
    ///
    /// assert_eq!(record.name(), b"r0");
    /// assert_eq!(record.sequence(), b"ATCG");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record).await
    }
}

async fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    record.clear();

    let mut len = match read_name(reader, record.name_mut()).await? {
        0 => return Ok(0),
        n => n,
    };

    len += read_line(reader, record.sequence_mut()).await?;
    len += consume_line(reader).await?;
    len += read_line(reader, record.quality_scores_mut()).await?;

    Ok(len)
}

async fn read_name<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    match reader.read_u8().await {
        Ok(b'@') => read_line(reader, buf).await.map(|n| n + 1),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid name prefix",
        )),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

async fn consume_line<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    let mut n = 0;
    let mut is_eol = false;

    while !is_eol {
        let buf = reader.fill_buf().await?;

        if buf.is_empty() {
            break;
        }

        let len = match buf.iter().position(|&b| b == LINE_FEED) {
            Some(i) => {
                is_eol = true;
                i + 1
            }
            None => buf.len(),
        };

        reader.consume(len);

        n += len;
    }

    Ok(n)
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

    #[tokio::test]
    async fn test_read_name() -> io::Result<()> {
        let data = b"@r0\n";
        let mut reader = &data[..];

        let mut buf = Vec::new();
        read_name(&mut reader, &mut buf).await?;

        assert_eq!(buf, b"r0");

        Ok(())
    }

    #[tokio::test]
    async fn test_consume_line() -> io::Result<()> {
        async fn t(mut data: &[u8], expected: &[u8]) -> io::Result<()> {
            consume_line(&mut data).await?;
            assert_eq!(data, expected);
            Ok(())
        }

        t(b"nd\nls\n", b"ls\n").await?;
        t(b"\nls\n", b"ls\n").await?;
        t(b"", b"").await?;

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
