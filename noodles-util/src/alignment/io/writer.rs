//! Alignment writer.

pub mod builder;
mod inner;

use std::io::{self, Write};

use noodles_sam::{self as sam, alignment::Record};

pub use self::builder::Builder;
use self::inner::Inner;

/// An alignment writer.
pub struct Writer<W>(Inner<W>)
where
    W: Write;

impl<W> Writer<W>
where
    W: Write,
{
    /// Writes a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::{self, io::Format};
    ///
    /// let mut writer = alignment::io::writer::Builder::default()
    ///     .set_format(Format::Bam)
    ///     .build_from_writer(io::sink())?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.0.write_header(header)
    }

    /// Writes an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    /// use noodles_util::alignment::{self, io::Format};
    ///
    /// let mut writer = alignment::io::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_writer(io::sink())?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = RecordBuf::default();
    /// writer.write_record(&header, &record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record<R>(&mut self, header: &sam::Header, record: &R) -> io::Result<()>
    where
        R: Record,
    {
        self.0.write_record(header, record)
    }

    /// Shuts down the alignment format writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::{self, io::Format};
    ///
    /// let mut writer = alignment::io::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_writer(io::sink())?;
    ///
    /// let header = sam::Header::default();
    /// writer.finish(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.0.finish(header)
    }
}
