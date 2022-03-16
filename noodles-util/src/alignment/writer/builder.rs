use std::io::Write;

use noodles_bam as bam;
use noodles_sam as sam;

use super::Writer;
use crate::alignment::Format;

/// An alignment writer builder.
pub struct Builder<W> {
    inner: W,
    format: Format,
}

impl<W> Builder<W>
where
    W: Write + 'static,
{
    pub(super) fn new(inner: W) -> Self {
        Self {
            inner,
            format: Format::Sam,
        }
    }

    /// Set the format of the output.
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Builds an alignment writer.
    pub fn build(self) -> Writer {
        let inner: Box<dyn sam::AlignmentWriter> = match self.format {
            Format::Sam => Box::new(sam::Writer::new(self.inner)),
            Format::Bam => Box::new(bam::Writer::new(self.inner)),
            Format::Cram => todo!(),
        };

        Writer { inner }
    }
}
