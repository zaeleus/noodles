mod builder;

pub use self::builder::Builder;

use std::io::{self, Write};

use noodles_sam as sam;

/// An alignment writer.
pub struct Writer {
    inner: Box<dyn sam::AlignmentWriter>,
}

impl Writer {
    /// Creates an alignment writer builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment;
    /// let builder = alignment::Writer::builder(Vec::new());
    /// ```
    pub fn builder<W>(inner: W) -> Builder<W>
    where
        W: Write + 'static,
    {
        Builder::new(inner)
    }

    /// Writes a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let mut writer = alignment::Writer::builder(Vec::new())
    ///     .set_format(Format::Sam)
    ///     .build();
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.inner.write_alignment_header(header)
    }

    /// Writes an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let mut writer = alignment::Writer::builder(Vec::new())
    ///     .set_format(Format::Sam)
    ///     .build();
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = sam::Record::default();
    /// writer.write_record(&header, &record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::AlignmentRecord,
    ) -> io::Result<()> {
        self.inner.write_alignment_record(header, record)
    }

    /// Shuts down the alignment format writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let mut writer = alignment::Writer::builder(Vec::new())
    ///     .set_format(Format::Sam)
    ///     .build();
    ///
    /// let header = sam::Header::default();
    /// writer.finish(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.inner.finish(header)
    }
}
