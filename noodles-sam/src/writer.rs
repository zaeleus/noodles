//! SAM writer.

mod builder;
mod num;
mod record;

use std::io::{self, Write};

pub use self::builder::Builder;
pub(crate) use self::record::write_record;
use super::{alignment::RecordBuf, Header};

/// A SAM writer.
///
/// The SAM format is comprised of two parts: 1) a header and 2) a list of records.
///
/// Each header line is prefixed with an `@` (at sign). The header is optional and may be empty.
///
/// SAM records are line-based and follow directly after the header or the start of the file until
/// EOF.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_sam::{self as sam, alignment::RecordBuf};
///
/// let mut writer = sam::Writer::new(Vec::new());
///
/// let header = sam::Header::builder().add_comment("noodles-sam").build();
/// writer.write_header(&header)?;
///
/// let record = RecordBuf::default();
/// writer.write_record(&header, &record)?;
///
/// let expected = b"@CO\tnoodles-sam
/// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
/// ";
///
/// assert_eq!(&writer.get_ref()[..], &expected[..]);
/// # Ok::<(), io::Error>(())
/// ```
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a SAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::Writer::new(Vec::new());
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
    /// use noodles_sam as sam;
    /// let mut writer = sam::Writer::new(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a SAM header.
    ///
    /// The SAM header is optional, though recommended to include. A call to this method can be
    /// omitted if it is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// let mut writer = sam::Writer::new(Vec::new());
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// writer.write_header(&header)?;
    /// assert_eq!(writer.get_ref(), b"@CO\tnoodles-sam\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write!(self.inner, "{header}")
    }

    /// Writes a SAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    ///
    /// let mut writer = sam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = RecordBuf::default();
    /// writer.write_record(&header, &record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &Header, record: &RecordBuf) -> io::Result<()> {
        write_record(&mut self.inner, header, record)
    }
}

impl<W> crate::alignment::io::Writer for Writer<W>
where
    W: Write,
{
    fn write_alignment_header(&mut self, header: &Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_alignment_record(
        &mut self,
        header: &Header,
        record: &dyn crate::alignment::Record,
    ) -> io::Result<()> {
        write_record(&mut self.inner, header, record)
    }

    fn finish(&mut self, _: &Header) -> io::Result<()> {
        Ok(())
    }
}
