use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

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

    /// Reads the raw SAM header.
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

    /// Reads a raw SAM record.
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
    /// reader.read_header().await?;
    ///
    /// let mut buf = String::new();
    /// reader.read_record(&mut buf).await?;
    /// assert_eq!(buf, "*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf).await
    }
}

async fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncBufRead + Unpin,
{
    const HEADER_PREFIX: u8 = b'@';
    const LINE_FEED: u8 = b'\n';

    let mut header_buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let buf = reader.fill_buf().await?;

        if (i == 0 || is_eol) && buf.first().map(|&b| b != HEADER_PREFIX).unwrap_or(true) {
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

async fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    let result = reader.read_line(buf).await;
    buf.pop();
    result
}
