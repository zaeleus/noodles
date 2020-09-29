use std::io::{self, Write};

use super::{record, Directive, Record};

/// A GFF writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a GFF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let writer = gff::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let writer = gff::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a GFF directive.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    ///
    /// let mut writer = gff::Writer::new(Vec::new());
    ///
    /// let version = gff::Directive::GffVersion(Default::default());
    /// writer.write_directive(&version)?;
    ///
    /// assert_eq!(writer.get_ref(), b"##gff-version 3\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_directive(&mut self, directive: &Directive) -> io::Result<()> {
        writeln!(self.inner, "{}", directive)
    }

    /// Writes a GFF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    ///
    /// let mut writer = gff::Writer::new(Vec::new());
    ///
    /// let version = gff::Directive::GffVersion(Default::default());
    /// writer.write_directive(&version)?;
    ///
    /// let record = gff::Record::default();
    /// writer.write_record(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write!(
            self.inner,
            "{seqid}\t{source}\t{ty}\t{start}\t{end}",
            seqid = record.reference_sequence_name(),
            source = record.source(),
            ty = record.ty(),
            start = record.start(),
            end = record.end(),
        )?;

        if let Some(score) = record.score() {
            write!(self.inner, "\t{}", score)?;
        } else {
            write!(self.inner, "\t{}", record::NULL_FIELD)?;
        }

        write!(self.inner, "\t{}", record.strand())?;

        if let Some(phase) = record.phase() {
            write!(self.inner, "\t{}", phase)?;
        } else {
            write!(self.inner, "\t{}", record::NULL_FIELD)?;
        }

        if record.attributes().is_empty() {
            writeln!(self.inner, "\t{}", record::NULL_FIELD)
        } else {
            writeln!(self.inner, "\t{}", record.attributes())
        }
    }
}
