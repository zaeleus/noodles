mod record;

use tokio::io::{self, AsyncBufRead};

use self::record::read_record;
use crate::fai::{Index, Record};

/// An async FASTA index reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Creates an async FASTA index reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let data = [];
    /// let mut reader = fai::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a FASTA index.
    ///
    /// The position of the stream is expected to be at the start or at the start of a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let data = b"sq0\t13\t5\t80\t81\nsq1\t21\t19\t80\t81\n";
    /// let mut reader = fai::r#async::io::Reader::new(&data[..]);
    ///
    /// let index = reader.read_index().await?;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// assert_eq!(index, fai::Index::from(vec![
    ///     fai::Record::new("sq0", 13, 5, line_base_count, line_width),
    ///     fai::Record::new("sq1", 21, 19, line_base_count, line_width),
    /// ]));
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        let mut buf = Vec::new();
        let mut records = Vec::new();

        loop {
            buf.clear();

            let mut record = Record::default();

            match read_record(&mut self.inner, &mut buf, &mut record).await? {
                0 => break,
                _ => records.push(record),
            }
        }

        Ok(Index::from(records))
    }
}
