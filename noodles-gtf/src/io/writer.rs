use std::io::{self, Write};

use crate::{Line, Record};

/// A GTF writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a GTF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let writer = gtf::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let writer = gtf::io::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Write};
    /// use noodles_gtf as gtf;
    /// let mut writer = gtf::io::Writer::new(Vec::new());
    /// writer.get_mut().write_all(b"ndls")?;
    /// assert_eq!(writer.get_ref(), b"ndls");
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let writer = gtf::io::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a [`Line`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// use gtf::line::Line;
    ///
    /// let mut writer = gtf::io::Writer::new(Vec::new());
    ///
    /// let version = Line::Comment(String::from("#format: gtf"));
    /// writer.write_line(&version)?;
    ///
    /// let comment = Line::Comment(String::from("noodles"));
    /// writer.write_line(&comment)?;
    ///
    /// let record = Line::Record(gtf::Record::default());
    /// writer.write_line(&record)?;
    ///
    /// let expected = b"##format: gtf
    /// #noodles
    /// .\t.\t.\t1\t1\t.\t.\t.\t
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    pub fn write_line(&mut self, line: &Line) -> io::Result<()> {
        writeln!(self.inner, "{line}")
    }

    /// Writes a GTF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    ///
    /// let mut writer = gtf::io::Writer::new(Vec::new());
    ///
    /// let record = gtf::Record::default();
    /// writer.write_record(&record)?;
    ///
    /// let expected = b".\t.\t.\t1\t1\t.\t.\t.\t\n";
    /// assert_eq!(writer.into_inner(), expected);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        writeln!(self.inner, "{record}")
    }
}
