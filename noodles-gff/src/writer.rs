use std::io::{self, Write};

use super::{Directive, Record};

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
        writeln!(self.inner, "{}", record)
    }
}
