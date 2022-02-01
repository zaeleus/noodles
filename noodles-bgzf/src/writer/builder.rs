use std::io::Write;

use super::Writer;
use crate::block;

/// A BGZF writer builder.
#[derive(Debug)]
pub struct Builder<W> {
    inner: W,
}

impl<W> Builder<W>
where
    W: Write,
{
    pub(super) fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Builds an async BGZF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::AsyncWriter::builder(Vec::new()).build();
    /// ```
    pub fn build(self) -> Writer<W> {
        Writer {
            inner: Some(self.inner),
            buf: Vec::with_capacity(block::MAX_UNCOMPRESSED_DATA_LENGTH),
        }
    }
}
