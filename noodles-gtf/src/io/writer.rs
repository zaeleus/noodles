mod line;

use std::io::{self, Write};

use self::line::write_line;
use crate::{LineBuf, RecordBuf};

/// A GTF writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// let writer = gtf::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// let mut writer = gtf::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf as gtf;
    /// let writer = gtf::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
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

    /// Writes a [`Line`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gtf::{self as gtf, LineBuf};
    ///
    /// let mut writer = gtf::io::Writer::new(Vec::new());
    ///
    /// let version = LineBuf::Comment(String::from("#format: gtf"));
    /// writer.write_line(&version)?;
    ///
    /// let comment = LineBuf::Comment(String::from("noodles"));
    /// writer.write_line(&comment)?;
    ///
    /// let record = LineBuf::Record(gtf::RecordBuf::default());
    /// writer.write_line(&record)?;
    ///
    /// let expected = b"##format: gtf
    /// #noodles
    /// .\t.\t.\t1\t1\t.\t.\t.\t
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    pub fn write_line(&mut self, line: &LineBuf) -> io::Result<()> {
        write_line(&mut self.inner, line)
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
    /// let record = gtf::RecordBuf::default();
    /// writer.write_record(&record)?;
    ///
    /// let expected = b".\t.\t.\t1\t1\t.\t.\t.\t\n";
    /// assert_eq!(writer.into_inner(), expected);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &RecordBuf) -> io::Result<()> {
        writeln!(self.inner, "{record}")
    }
}
