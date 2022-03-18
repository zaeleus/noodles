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
    pub fn builder<W>(inner: W) -> Builder<W>
    where
        W: Write + 'static,
    {
        Builder::new(inner)
    }

    /// Writes a SAM header.
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.inner.write_alignment_header(header)
    }

    /// Writes an alignment record.
    pub fn write_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::AlignmentRecord,
    ) -> io::Result<()> {
        self.inner.write_alignment_record(header, record)
    }

    /// Shuts down the alignment format writer.
    pub fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.inner.finish(header)
    }
}
