mod line;

use std::{
    io::{self, BufRead, Read, Seek},
    iter,
};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::{self as csi, BinningIndex};
use noodles_gff::feature::RecordBuf;

use crate::{Line, LineBuf, Record};

/// A GTF reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// let reader = gtf::io::Reader::new(io::empty());
    /// let _ = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// let mut reader = gtf::io::Reader::new(io::empty());
    /// let _ = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// let reader = gtf::io::Reader::new(io::empty());
    /// let _ = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
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
    /// let reader = gtf::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a GTF line.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    ///
    /// let src = b"##format: gtf\nsq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tid 0;\n";
    /// let mut reader = gtf::io::Reader::new(&src[..]);
    ///
    /// let mut line = gtf::Line::default();
    ///
    /// reader.read_line(&mut line)?;
    /// assert_eq!(line.kind(), gtf::line::Kind::Comment);
    ///
    /// reader.read_line(&mut line)?;
    /// assert_eq!(line.kind(), gtf::line::Kind::Record);
    ///
    /// assert_eq!(reader.read_line(&mut line)?, 0);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_line(&mut self, line: &mut Line) -> io::Result<usize> {
        line::read_line(&mut self.inner, line)
    }

    /// Returns an iterator over lines starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    ///
    /// let src = b"##format: gtf\nsq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tid 0;\n";
    /// let mut reader = gtf::io::Reader::new(&src[..]);
    ///
    /// let mut lines = reader.lines();
    ///
    /// let line = lines.next().transpose()?;
    /// assert_eq!(line.map(|l| l.kind()), Some(gtf::line::Kind::Comment));
    ///
    /// let line = lines.next().transpose()?;
    /// assert_eq!(line.map(|l| l.kind()), Some(gtf::line::Kind::Record));
    ///
    /// assert!(lines.next().is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn lines(&mut self) -> impl Iterator<Item = io::Result<Line>> + '_ {
        let mut line = Line::default();

        iter::from_fn(move || match self.read_line(&mut line) {
            Ok(0) => None,
            Ok(_) => Some(Ok(line.clone())),
            Err(e) => Some(Err(e)),
        })
    }

    /// Returns an iterator over line buffers starting from the current stream position.
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
    /// let mut reader = gtf::io::Reader::new(&data[..]);
    ///
    /// let mut lines = reader.line_bufs();
    ///
    /// let line = lines.next().transpose()?;
    /// assert_eq!(line, Some(gtf::LineBuf::Comment(String::from("#format: gtf"))));
    ///
    /// let line = lines.next().transpose()?;
    /// assert!(matches!(line, Some(gtf::LineBuf::Record(_))));
    ///
    /// assert!(lines.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn line_bufs(&mut self) -> impl Iterator<Item = io::Result<LineBuf>> + '_ {
        let mut line = Line::default();

        iter::from_fn(move || match self.read_line(&mut line) {
            Ok(0) => None,
            Ok(_) => Some(
                LineBuf::try_from(line.clone())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            ),
            Err(e) => Some(Err(e)),
        })
    }

    /// Returns an iterator over record buffers starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    ///
    /// let src = b"##format: gtf
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"g0\"; transcript_id \"t0\";
    /// ";
    /// let mut reader = gtf::io::Reader::new(&src[..]);
    ///
    /// let mut record_bufs = reader.record_bufs();
    ///
    /// let record = record_bufs.next().transpose()?;
    /// assert_eq!(record.map(|r| r.start()), Position::new(8));
    /// // ...
    ///
    /// assert!(record_bufs.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn record_bufs(&mut self) -> impl Iterator<Item = io::Result<RecordBuf>> + '_ {
        let mut lines = self.line_bufs();

        iter::from_fn(move || loop {
            match lines.next()? {
                Ok(LineBuf::Record(r)) => return Some(Ok(r)),
                Ok(_) => {}
                Err(e) => return Some(Err(e)),
            }
        })
    }
}

impl<R> Reader<bgzf::io::Reader<R>>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'r, I>(
        &'r mut self,
        index: &I,
        region: &'r Region,
    ) -> io::Result<impl Iterator<Item = io::Result<RecordBuf>> + 'r>
    where
        I: BinningIndex,
    {
        let header = index
            .header()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index header"))?;

        let reference_sequence_id = header
            .reference_sequence_names()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "missing reference sequence name",
                )
            })?;

        let chunks = index.query(reference_sequence_id, region.interval())?;

        let records = csi::io::Query::new(&mut self.inner, chunks)
            .indexed_records(header)
            .filter_by_region(region)
            .map(|result| {
                result.and_then(|r| {
                    Record::try_new(r.as_ref())
                        .and_then(|record| RecordBuf::try_from_feature_record(&record))
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
            });

        Ok(records)
    }
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    match reader.read_line(buf)? {
        0 => Ok(0),
        n => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
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
