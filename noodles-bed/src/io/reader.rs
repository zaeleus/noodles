use std::{
    io::{self, BufRead},
    iter,
    str::FromStr,
};

use crate::feature::RecordBuf;

/// A BED reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::new(&data[..]);
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
    /// use noodles_bed as bed;
    /// let data = [];
    /// let mut reader = bed::io::Reader::new(&data[..]);
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
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a BED reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let data = [];
    /// let reader = bed::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a raw BED line.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    ///
    /// let data = b"sq0\t8\t13\n";
    /// let mut reader = bed::io::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_line(&mut buf)?;
    ///
    /// assert_eq!(buf, "sq0\t8\t13");
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let data = b"sq0\t7\t13\n# sq0\t20\t34\n";
    /// let mut reader = bed::io::Reader::new(&data[..]);
    ///
    /// let mut records = reader.records::<3>();
    ///
    /// let record = records.next().transpose()?;
    /// assert_eq!(record.map(|r| r.start_position()), Position::new(8));
    /// // ...
    ///
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<const N: u8>(&mut self) -> impl Iterator<Item = io::Result<RecordBuf<N>>> + '_
    where
        RecordBuf<N>: FromStr<Err = crate::feature::record_buf::ParseError>,
    {
        const COMMENT_PREFIX: char = '#';

        let mut buf = String::new();

        iter::from_fn(move || loop {
            buf.clear();

            match self.read_line(&mut buf) {
                Ok(0) => return None,
                Ok(_) => {
                    if buf.starts_with(COMMENT_PREFIX) {
                        continue;
                    } else {
                        return Some(
                            buf.parse()
                                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                        );
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        })
    }
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    match reader.read_line(buf) {
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

    #[test]
    fn test_read_line() -> io::Result<()> {
        fn t(buf: &mut String, mut reader: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut reader, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = String::new();

        t(&mut buf, b"noodles\n", "noodles")?;
        t(&mut buf, b"noodles\r\n", "noodles")?;
        t(&mut buf, b"noodles", "noodles")?;

        Ok(())
    }
}
