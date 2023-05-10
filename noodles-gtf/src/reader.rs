use std::{
    io::{self, BufRead},
    iter,
};

use super::{Line, Record};

/// A GTF reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a GTF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let data = [];
    /// let reader = gtf::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a raw GTF line.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    ///
    /// let data = b"sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"ndls0\"; transcript_id \"ndls0\";";
    /// let mut reader = gtf::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_line(&mut buf)?;
    ///
    /// assert_eq!(
    ///     buf,
    ///     "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"ndls0\"; transcript_id \"ndls0\";"
    /// );
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf)
    }

    /// Returns an iterator over lines starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    ///
    /// let data = b"##format: gtf
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"g0\"; transcript_id \"t0\";
    /// ";
    /// let mut reader = gtf::Reader::new(&data[..]);
    ///
    /// let mut lines = reader.lines();
    ///
    /// let line = lines.next().transpose()?;
    /// assert_eq!(line, Some(gtf::Line::Comment(String::from("#format: gtf"))));
    ///
    /// let line = lines.next().transpose()?;
    /// assert!(matches!(line, Some(gtf::Line::Record(_))));
    ///
    /// assert!(lines.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn lines(&mut self) -> impl Iterator<Item = io::Result<Line>> + '_ {
        let mut buf = String::new();

        iter::from_fn(move || {
            buf.clear();

            match self.read_line(&mut buf) {
                Ok(0) => None,
                Ok(_) => Some(
                    buf.parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                ),
                Err(e) => Some(Err(e)),
            }
        })
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    ///
    /// let data = b"##format: gtf
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"g0\"; transcript_id \"t0\";
    /// ";
    /// let mut reader = gtf::Reader::new(&data[..]);
    ///
    /// let mut records = reader.records();
    ///
    /// let record = records.next().transpose()?;
    /// assert_eq!(record.map(|r| r.start()), Position::new(8));
    /// // ...
    ///
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records(&mut self) -> impl Iterator<Item = io::Result<Record>> + '_ {
        let mut lines = self.lines();

        iter::from_fn(move || loop {
            match lines.next()? {
                Ok(Line::Record(r)) => return Some(Ok(r)),
                Ok(_) => {}
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
        fn t(buf: &mut String, mut src: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut src, buf)?;
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
