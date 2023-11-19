//! GFF reader and iterators.

mod lines;
mod records;

pub use self::{lines::Lines, records::Records};

use std::{
    io::{self, BufRead, Read, Seek},
    mem, str,
};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::{self as csi, BinningIndex};

use super::{lazy, Record};

const LINE_FEED: char = '\n';
const CARRIAGE_RETURN: char = '\r';

/// A GFF reader.
pub struct Reader<R> {
    inner: R,
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
    /// use noodles_gff as gff;
    /// let data = b"##gff-version 3\n";
    /// let mut reader = gff::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    ///
    /// let data = b"##gff-version 3\n";
    /// let reader = gff::Reader::new(&data[..]);
    ///
    /// let _ = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    ///
    /// let data = b"##gff-version 3
    /// #format: gff3
    /// ";
    /// let mut reader = gff::Reader::new(&data[..]);
    /// reader.read_line(&mut String::new())?;
    ///
    /// assert_eq!(reader.into_inner(), b"#format: gff3\n");
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads a raw GFF line.
    ///
    /// This reads from the underlying stream until a newline is reached and appends it to the
    /// given buffer, sans the final newline character. The buffer can subsequently be parsed as a
    /// [`crate::Line`].
    ///
    /// It is more ergonomic to read records using an iterator (see [`Self::lines`]), but using
    /// this method allows control of the line buffer and whether the raw line should be parsed.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
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
    /// let mut reader = gff::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_line(&mut buf)?;
    /// assert_eq!(buf, "##gff-version 3");
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf)
    }

    /// Returns an iterator over lines starting from the current stream position.
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
    /// use noodles_gff as gff;
    ///
    /// let data = b"##gff-version 3
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0
    /// ";
    /// let mut reader = gff::Reader::new(&data[..]);
    /// let mut lines = reader.lines();
    ///
    /// let line = lines.next().transpose()?;
    /// assert!(matches!(line, Some(gff::Line::Directive(_))));
    ///
    /// let line = lines.next().transpose()?;
    /// assert!(matches!(line, Some(gff::Line::Record(_))));
    ///
    /// assert!(lines.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn lines(&mut self) -> Lines<'_, R> {
        Lines::new(self)
    }

    /// Reads a single line without eagerly decoding it.
    pub fn read_lazy_line(&mut self, line: &mut lazy::Line) -> io::Result<usize> {
        const DEFAULT_LINE: lazy::Line = lazy::Line::Comment(String::new());
        const DIRECTIVE_PREFIX: &str = "##";

        let prev_line = mem::replace(line, DEFAULT_LINE);
        let mut buf = prev_line.into();

        match peek_line_type(&mut self.inner)? {
            Some(LineType::Comment) => {
                let n = read_line(&mut self.inner, &mut buf)?;

                *line = if buf.starts_with(DIRECTIVE_PREFIX) {
                    lazy::Line::Directive(buf)
                } else {
                    lazy::Line::Comment(buf)
                };

                Ok(n)
            }
            Some(LineType::Record) => {
                let (n, bounds) = read_lazy_record(&mut self.inner, &mut buf)?;
                *line = lazy::Line::Record(lazy::Record { buf, bounds });
                Ok(n)
            }
            None => Ok(0),
        }
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
    /// let mut reader = gff::Reader::new(&data[..]);
    /// let mut records = reader.records();
    ///
    /// assert!(records.next().transpose()?.is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self.lines())
    }
}

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'r, I>(
        &'r mut self,
        index: &I,
        region: &'r Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + 'r>
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
                    r.as_ref()
                        .parse()
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

enum LineType {
    Comment,
    Record,
}

fn peek_line_type<R>(reader: &mut R) -> io::Result<Option<LineType>>
where
    R: BufRead,
{
    const COMMENT_PREFIX: u8 = b'#';

    let src = reader.fill_buf()?;

    Ok(src.first().map(|&b| match b {
        COMMENT_PREFIX => LineType::Comment,
        _ => LineType::Record,
    }))
}

fn read_lazy_record<R>(
    reader: &mut R,
    buf: &mut String,
) -> io::Result<(usize, lazy::record::Bounds)>
where
    R: BufRead,
{
    buf.clear();

    let mut len = 0;
    let mut bounds = lazy::record::Bounds::default();

    len += read_field(reader, buf)?;
    bounds.reference_sequence_name_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.source_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.type_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.start_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.end_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.score_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.strand_end = buf.len();

    len += read_field(reader, buf)?;
    bounds.phase_end = buf.len();

    len += read_line(reader, buf)?;

    Ok((len, bounds))
}

fn read_field<R>(reader: &mut R, dst: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const DELIMITER: u8 = b'\t';

    let mut is_delimiter = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if is_delimiter || src.is_empty() {
            break;
        }

        let (buf, n) = match src.iter().position(|&b| b == DELIMITER) {
            Some(i) => {
                is_delimiter = true;
                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        let s = str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        dst.push_str(s);

        len += n;

        reader.consume(n);
    }

    Ok(len)
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

        for result in reader.records() {
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

        for result in reader.records() {
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
