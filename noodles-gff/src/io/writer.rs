mod line;

use std::io::{self, Write};

use self::line::write_line;
use crate::{feature::RecordBuf, DirectiveBuf, LineBuf};

/// A GFF writer.
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
    /// use noodles_gff as gff;
    /// let writer = gff::io::Writer::new(io::sink());
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
    /// use noodles_gff as gff;
    /// let mut writer = gff::io::Writer::new(io::sink());
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
    /// use noodles_gff as gff;
    /// let writer = gff::io::Writer::new(io::sink());
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
    /// Creates a GFF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let writer = gff::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a [`LineBuf`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use bstr::BString;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}, LineBuf};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = LineBuf::Directive(gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// ));
    /// writer.write_line(&version)?;
    ///
    /// let comment = LineBuf::Comment(BString::from("noodles"));
    /// writer.write_line(&comment)?;
    ///
    /// let record = LineBuf::Record(gff::feature::RecordBuf::default());
    /// writer.write_line(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// #noodles
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    pub fn write_line(&mut self, line: &LineBuf) -> io::Result<()> {
        write_line(&mut self.inner, line)
    }

    /// Writes a GFF directive.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// );
    /// writer.write_directive(&version)?;
    ///
    /// assert_eq!(writer.get_ref(), b"##gff-version 3\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_directive(&mut self, directive: &DirectiveBuf) -> io::Result<()> {
        line::write_directive(&mut self.inner, directive)?;
        line::write_newline(&mut self.inner)
    }

    /// Writes a GFF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// );
    /// writer.write_directive(&version)?;
    ///
    /// let record = gff::feature::RecordBuf::default();
    /// writer.write_record(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &RecordBuf) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// );
    /// writer.write_directive(&version)?;
    ///
    /// let record = gff::feature::RecordBuf::default();
    /// writer.write_feature_record(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_feature_record(&mut self, record: &dyn crate::feature::Record) -> io::Result<()> {
        line::write_record(&mut self.inner, record)?;
        line::write_newline(&mut self.inner)
    }
}
