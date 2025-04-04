//! GFF reader and iterators.

pub(crate) mod line;
mod line_bufs;
mod lines;
mod record_bufs;

pub use self::{line_bufs::LineBufs, lines::Lines, record_bufs::RecordBufs};

use std::io::{self, BufRead, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::{self as csi, BinningIndex};

use crate::{feature::RecordBuf, Line};

/// A GFF reader.
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
    /// use noodles_gff as gff;
    /// let reader = gff::io::Reader::new(io::empty());
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
    /// use noodles_gff as gff;
    /// let mut reader = gff::io::Reader::new(io::empty());
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
    /// use noodles_gff as gff;
    /// let reader = gff::io::Reader::new(io::empty());
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
    /// Creates a GFF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    /// let reader = gff::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns an iterator over line buffers starting from the current stream position.
    ///
    /// When using this, the caller is responsible to stop reading at either EOF or when the
    /// `FASTA` directive is read, whichever comes first.
    ///
    /// Unlike [`Self::read_line`], each line is parsed as a [`crate::Line`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, LineBuf};
    ///
    /// let data = b"##gff-version 3
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0
    /// ";
    /// let mut reader = gff::io::Reader::new(&data[..]);
    /// let mut lines = reader.line_bufs();
    ///
    /// let line = lines.next().transpose()?;
    /// assert!(matches!(line, Some(LineBuf::Directive(_))));
    ///
    /// let line = lines.next().transpose()?;
    /// assert!(matches!(line, Some(LineBuf::Record(_))));
    ///
    /// assert!(lines.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn line_bufs(&mut self) -> LineBufs<'_, R> {
        LineBufs::new(self)
    }

    /// Reads a single line without eagerly decoding it.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    ///
    /// let data = b"##gff-version 3\n";
    /// let mut reader = gff::io::Reader::new(&data[..]);
    ///
    /// let mut line = gff::Line::default();
    ///
    /// reader.read_line(&mut line)?;
    /// assert_eq!(line.kind(), gff::line::Kind::Directive);
    ///
    /// assert_eq!(reader.read_line(&mut line)?, 0);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_line(&mut self, line: &mut Line) -> io::Result<usize> {
        line::read_line(&mut self.inner, line)
    }

    /// Returns an iterator over lines starting from the current stream position.
    ///
    /// When using this, the caller is responsible to stop reading at either EOF or when the
    /// `FASTA` directive is read, whichever comes first.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::key};
    ///
    /// let mut reader = gff::io::Reader::new(io::empty());
    ///
    /// for result in reader.lines() {
    ///     let line = result?;
    ///
    ///     if let Some(key::FASTA) = line.as_directive().map(|directive| directive.key().as_ref()) {
    ///         break;
    ///     }
    ///
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn lines(&mut self) -> Lines<'_, R> {
        Lines::new(self)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// This filters lines for only records. It stops at either EOF or when the `FASTA` directive
    /// is read, whichever comes first.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    ///
    /// let data = b"##gff-version 3
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0
    /// ";
    /// let mut reader = gff::io::Reader::new(&data[..]);
    /// let mut records = reader.record_bufs();
    ///
    /// assert!(records.next().transpose()?.is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn record_bufs(&mut self) -> RecordBufs<'_, R> {
        RecordBufs::new(self.line_bufs())
    }
}

impl<R> Reader<bgzf::io::Reader<R>>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi as csi;
    /// use noodles_gff as gff;
    ///
    /// let mut reader = File::open("annotations.gff3.gz")
    ///     .map(bgzf::io::Reader::new)
    ///     .map(gff::io::Reader::new)?;
    ///
    /// let index = csi::fs::read("annotations.gff3.gz.csi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
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
                    let line = Line(r.as_ref().into());

                    line.as_record()
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidData, "line is not a record")
                        })?
                        .and_then(|record| RecordBuf::try_from_feature_record(&record))
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
    fn test_records() -> io::Result<()> {
        let data = b"\
##gff-version 3
sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0
";

        let mut reader = Reader::new(&data[..]);
        let mut n = 0;

        for result in reader.record_bufs() {
            let _ = result?;
            n += 1;
        }

        assert_eq!(n, 1);

        Ok(())
    }

    #[test]
    fn test_records_with_fasta_directive() -> io::Result<()> {
        let data = b"\
##gff-version 3
sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0
##FASTA
>sq0
ACGT
";

        let mut reader = Reader::new(&data[..]);
        let mut n = 0;

        for result in reader.record_bufs() {
            let _ = result?;
            n += 1;
        }

        assert_eq!(n, 1);

        Ok(())
    }

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
