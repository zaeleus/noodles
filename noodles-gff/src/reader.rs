mod lines;
mod records;

pub use self::{lines::Lines, records::Records};

use std::io::{self, BufRead};

pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
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
    /// [`gff::Line`].
    ///
    /// It is more ergonomic to read records using an iterator (see [`lines`]), but using this
    /// method allows control of the line buffer and whether the raw line should be parsed.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// [`gff::Line`]: struct.Line.html
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
        let result = self.inner.read_line(buf);
        buf.pop();
        result
    }

    /// Returns an iterator over lines starting from the current stream position.
    ///
    /// This stops at either EOF or when the `FASTA` directive is read, whichever comes first.
    ///
    /// Unlike [`read_line`], each line is parsed as a [`gff::Line`].
    ///
    /// [`read_line`]: #method.read_line
    /// [`gff::Line`]: struct.Line.html
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
    pub fn lines(&mut self) -> Lines<R> {
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
    /// let mut reader = gff::Reader::new(&data[..]);
    /// let mut records = reader.records();
    ///
    /// assert!(records.next().transpose()?.is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<R> {
        Records::new(self.lines())
    }
}

#[cfg(test)]
mod tests {
    use crate::Line;

    use super::*;

    #[test]
    fn test_lines_with_fasta_directive() -> io::Result<()> {
        let data = b"\
##gff-version 3
sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0
##FASTA
>sq0
NNNNNNNNNNNNNNNN
";

        let mut reader = Reader::new(&data[..]);
        let mut lines = reader.lines();

        let line = lines.next().transpose()?;
        assert!(matches!(line, Some(Line::Directive(_))));

        let line = lines.next().transpose()?;
        assert!(matches!(line, Some(Line::Record(_))));

        assert!(lines.next().is_none());

        Ok(())
    }
}
